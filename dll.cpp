#include "tigr.h"
#include "math_101.h"
#include "draw.h"
#include "dll.h"
#include "routine.h"

#include <stdio.h>
#include <assert.h>

malloc_fn* g_malloc;
free_fn* g_free;
#define ALLOC g_malloc
#define FREE g_free

inline int key_pressed(int key) { return tigrKeyDown(screen, key); }
inline int key_down(int key) { return tigrKeyHeld(screen, key); }

float g_seconds = 0;
float g_prev_seconds = 0;
float g_dt = 0;

inline void update_time(float dt)
{
	g_prev_seconds = g_seconds;
	g_seconds += dt;
	g_dt = dt;
}

inline float delta_time()
{
	return g_dt;
}

inline bool on_interval(float interval)
{
	int prev = (int)(g_prev_seconds / interval);
	int next = (int)(g_seconds / interval);
	return prev < next;
}

inline bool between_interval(float interval)
{
	return fmodf(g_seconds, interval * 2) >= interval;
}

#define ROCKETS_MAX 16

struct RocketBarn
{
	int count;
	int capacity;
	bool alive[ROCKETS_MAX];
	float elapsed[ROCKETS_MAX];
	float hang_time[ROCKETS_MAX];
	Routine rt[ROCKETS_MAX];
	v2 p[ROCKETS_MAX];
	v2 c0[ROCKETS_MAX];
	v2 c1[ROCKETS_MAX];
	v2 target[ROCKETS_MAX];

	float prev_dt;
	int hit_count;
	v2 hits[ROCKETS_MAX];

	bool add(v2 start, v2 end, float duration);
	void update();
	void draw();
	bool try_pop_hit(v2* hit_out);
};

#define BULLETS_MAX 16

struct BulletBarn
{
	bool alive[BULLETS_MAX];
	v2 p[BULLETS_MAX];

	bool add();
	void update();
	void draw();
};

#define ASTEROIDS_MAX 256

struct AsteroidBarn
{
	bool alive[ASTEROIDS_MAX];
	float angular_velocity[ASTEROIDS_MAX];
	v2 velocity[ASTEROIDS_MAX];
	polygon poly[ASTEROIDS_MAX];
	v2 center_of_mass[ASTEROIDS_MAX];
	float slice_timeout[ASTEROIDS_MAX];

	void add(v2 p, v2 v);
	void add(polygon p, v2 v, float a, float timeout);
	void update();
	void slice(ray r);
	void draw();
};

#define TRAIL_MAX 128

struct TrailBarn
{
	bool alive[TRAIL_MAX];
	float lifespan[TRAIL_MAX];
	float time_left[TRAIL_MAX];
	aabb shape[TRAIL_MAX];

	void add(aabb box, float duration);
	void update();
	void draw();
};

struct SpriteDef
{
	const char* name;
	uint64_t hash;
	Routine co;
	int image_count;
	void (*load)(Tigr**);
	float durations[64];
	inline void play(const char* state_name) { nav_next(co, state_name); }
};

enum AnimationDirection
{
	ANIMATION_DIRECTION_FORWARD,
	ANIMATION_DIRECTION_BACKWARD,
	ANIMATION_DIRECTION_PINGPONG,
};

struct Sprite
{
	Sprite() { memset(this, 0, sizeof(*this)); animate_forwards = true; }
	Sprite(const char* name, AnimationDirection dir = ANIMATION_DIRECTION_FORWARD) { set(name, dir); }

	SpriteDef* def;
	v2 p;
	float time;
	bool paused;
	bool loop;
	bool animate_forwards;
	bool animate_backwards;
	bool animate_pingpong;
	bool finished;
	int frame;
	Tigr** images;

	void set(const char* name, AnimationDirection dir = ANIMATION_DIRECTION_FORWARD);
	bool valid();
	void update();
	void draw(TPixel tint = tigrRGBA(0xFF,0xFF,0xFF,0xFF));
	bool on_finished();
};

struct PlayerShip
{
	v2 p;
	aabb bounds;
	Routine rt_movement;
	Routine rt_weapons;
	Routine rt_trail;
	bool charging_laser;
	int charging_shot;
	bool fired_laser;
	bool firing_rockets;
	bool shielding;
	float shield_time;
	RocketBarn rockets;
	BulletBarn bullets;
	ray laser_trail;
	Sprite sprite;
	int facing_index;
	float booster_time;
	float turn_time;
	v2 dir;
	v2 old_dir;

	void reset();
};

#define ANIMS_MAX 256

struct AnimBarn
{
	bool alive[ANIMS_MAX];
	bool fade[ANIMS_MAX];
	float alpha[ANIMS_MAX];
	TPixel tint[ANIMS_MAX];
	Sprite sprite[ANIMS_MAX];

	inline void add(Sprite sprite, bool fade = false, TPixel tint = tigrRGBA(0xFF,0xFF,0xFF,0xFF))
	{
		for (int i = 0; i < ANIMS_MAX; ++i) {
			if (alive[i]) continue;
			alive[i] = true;
			this->fade[i] = fade;
			this->alpha[i] = 1.0f;
			this->tint[i] = tint;
			this->sprite[i] = sprite;
			break;
		}
	}

	inline void update()
	{
		for (int i = 0; i < ANIMS_MAX; ++i) {
			if (!alive[i]) continue;
			sprite[i].update();
			if (sprite[i].on_finished()) {
				alive[i] = false;
				continue;
			}
		}
	}

	inline void draw()
	{
		for (int i = 0; i < ANIMS_MAX; ++i) {
			if (!alive[i]) continue;
			sprite[i].draw(tint[i]);
		}
	}
};

#define ANIMATIONS_MAX 64
#define ANIMATION_FRAMES_MAX 16

struct game
{
	rnd_t rnd;
	PlayerShip player;
	AsteroidBarn asteroids;
	TrailBarn trails;
	AnimBarn anims;
	Tigr* images[ANIMATIONS_MAX][ANIMATION_FRAMES_MAX];
};

game* g;

// name - Name of the sprite.
// count - Number of images in the animations.
// ... - Array of floats for frame durations in seconds.
#define REGISTER_SPRITE(name, count, ...)                  \
	{                                                      \
		name,                                              \
		fnv1a(name),                                       \
		Routine(),                                         \
		count,                                             \
		[](Tigr** images) {                                \
			char buf[256];                                 \
			static_assert(count <= ANIMATION_FRAMES_MAX, "When using REGISTER_SPRITE count must be <= ANIMATION_FRAMES_MAX."); \
			for (int i = 0; i < count; ++i) {              \
				sprintf(buf, "art/%s%d.png", name, i + 1); \
				images[i] = tigrLoadImage(buf);            \
			}                                              \
		},                                                 \
		__VA_ARGS__, /* The frame durations. */            \
	}                                                      \

SpriteDef sprite_defs[] = {
	REGISTER_SPRITE(
		"explosion",
		5,
		{ 0.025f, 0.025f, 0.025f, 0.05f, 0.1f }
	),
	REGISTER_SPRITE(
		"charge",
		2,
		{ 0.1f }
	),
	REGISTER_SPRITE(
		"ship",
		14
	),
};

void load_all_sprites()
{
	for (int i = 0; i < sizeof(sprite_defs) / sizeof(*sprite_defs); ++i) {
		sprite_defs[i].load(g->images[i]);
	}
}

inline void Sprite::set(const char* name, AnimationDirection dir)
{
	*this = Sprite();
	if (dir == ANIMATION_DIRECTION_BACKWARD) {
		animate_forwards = false;
		animate_backwards = true;
	} else if (dir == ANIMATION_DIRECTION_PINGPONG) {
		animate_pingpong = true;
	}
	uint64_t h = fnv1a(name);
	for (int i = 0; i < sizeof(sprite_defs) / sizeof(*sprite_defs); ++i) {
		if (h == sprite_defs[i].hash) {
			def = sprite_defs + i;
			images = g->images[i];
		}
	}
}

inline bool Sprite::valid()
{
	if (!def) return false;
	if (frame == -1) return false;
	return true;
}

inline void Sprite::update()
{
	if (!valid()) return;
	if (paused) return;
	float dt = delta_time();
	time += dt;
	finished = false;
	if (time > def->durations[frame]) {
		if (animate_forwards) {
			assert(animate_backwards == false);
			if (frame + 1 == def->image_count) {
				if (animate_pingpong) {
					animate_forwards = false;
					animate_backwards = true;
					frame = max(0, frame - 1);
				} else {
					if (loop) frame = 0;
					finished = true;
				}
			} else {
				frame++;
			}
		} else if (animate_backwards) {
			if (frame - 1 < 0) {
				finished = true;
				if (animate_pingpong) {
					if (loop) {
						animate_forwards = true;
						animate_backwards = false;
						frame = min(def->image_count - 1, frame + 1);
					}
				} else {
					if (loop) frame = def->image_count - 1;
				}
			} else {
				frame--;
			}
		}
		if (loop && !finished) time = 0;
	}
}

inline void Sprite::draw(TPixel tint)
{
	if (!valid()) return;
	v2 pos = world_to_screen(p);
	int w = (int)(0.5f + images[frame]->w / 2.0f);
	int h = (int)(0.5f + images[frame]->h / 2.0f);
	tigrBlitTint(screen, images[frame], (int)(pos.x - w), (int)(pos.y - h), 0, 0, (int)(w*2), (int)(h*2), tint);
}

inline bool Sprite::on_finished()
{
	return finished;
}

v2 rnd_next_range(v2 lo, v2 hi) { return v2(rnd_next_range(g->rnd, lo.x, hi.x), rnd_next_range(g->rnd, lo.y, hi.y)); }

inline float fade(float t, float lo, float hi)
{
	return map(1.0f - smoothstep(t), lo, hi);
}

inline int fade_alpha(float t, float lo, float hi)
{
	return (int)(0.5f + map(1.0f - smoothstep(t), lo, hi) * 255.0f);
}

void TrailBarn::add(aabb box, float duration)
{
	for (int i = 0; i < TRAIL_MAX; ++i) {
		if (alive[i]) continue;
		alive[i] = true;
		shape[i] = box;
		lifespan[i] = duration;
		time_left[i] = duration;
		break;
	}
}

void TrailBarn::update()
{
	for (int i = 0; i < TRAIL_MAX; ++i) {
		if (!alive[i]) continue;
		time_left[i] -= delta_time();
		if (time_left[i] < 0) alive[i] = false;
	}
}

void TrailBarn::draw()
{
	for (int i = 0; i < TRAIL_MAX; ++i) {
		if (!alive[i]) continue;
		draw_box(shape[i], color_white(fade_alpha(1.0f - (time_left[i] / lifespan[i]), 0.4f, 1.0f)));
	}
}

const float player_shield_max = 3.0f;

void PlayerShip::reset()
{
	memset(this, 0, sizeof(*this));
	rockets.capacity = 5;
	sprite = Sprite("ship");
	facing_index = 3;
}

bool RocketBarn::add(v2 start, v2 end, float duration)
{
	if (count >= capacity) return false;
	for (int i = 0; i < ROCKETS_MAX; ++i) {
		if (!alive[i]) {
			alive[i] = true;
			hang_time[i] = duration;
			elapsed[i] = 0;
			rt[i] = Routine();
			v2 d = norm(end - start);
			p[i] = start;
			c0[i] = (start - d * 150.0f + skew(d) * rnd_next_range(g->rnd, 25.0f, 250.0f) * sign(rnd_next_range(g->rnd, -1.0f, 1.0f)));
			c1[i] = rnd_next_range((start + end) * 0.5f - v2(50,50), (start + end) * 0.5f+ v2(50,50));
			target[i] = end;
			++count;
			return true;
		}
	}
	return false;
}

void RocketBarn::update()
{
	prev_dt = delta_time();
	hit_count = 0;
	for (int i = 0; i < ROCKETS_MAX; ++i) {
		if (!alive[i]) continue;
		rt_begin(rt[i], delta_time());
		rt_seconds(hang_time[i])
		{
			elapsed[i] = rt.elapsed;
			if (on_interval(0.015f)) {
				float t = elapsed[i];
				v2 at = bezier(p[i], c0[i], c1[i], target[i], ease_in_sin(t));
				g->trails.add(aabb(at - v2(1,1), at + v2(1,1)), 0.15f);
			}
		}
		rt_once()
		{
			alive[i] = false;
			Sprite s = Sprite("explosion");
			s.p = target[i];
			g->anims.add(s);
			hits[hit_count++] = target[i];
		}
		rt_end();
	}
}

void RocketBarn::draw()
{
	for (int i = 0; i < ROCKETS_MAX; ++i) {
		if (!alive[i]) continue;
		float t0 = (elapsed[i] - prev_dt * 2 / hang_time[i]);
		float t1 = (elapsed[i] - prev_dt / hang_time[i]);
		float t = elapsed[i];
		v2 at0 = bezier(p[i], c0[i], c1[i], target[i], ease_in_sin(t0));
		v2 at1 = bezier(p[i], c0[i], c1[i], target[i], ease_in_sin(t1));
		v2 at = bezier(p[i], c0[i], c1[i], target[i], ease_in_sin(t));
		draw_circle_fill(circle(at, 3), color_white());
		draw_circle_fill(circle(at1, 4), color_white(0xBB));
		draw_circle_fill(circle(at0, 4), color_white(0x88));
	}
}

bool RocketBarn::try_pop_hit(v2* hit_out)
{
	if (hit_count) {
		*hit_out = hits[--hit_count];
		return true;
	} else {
		return false;
	}
}

bool BulletBarn::add()
{
	for (int i = 0; i < BULLETS_MAX; ++i) {
		if (!alive[i]) {
			alive[i] = true;
			p[i] = top(g->player.bounds) + v2(0, 8);
			return true;
		}
	}
	return false;
}

void BulletBarn::update()
{
	for (int i = 0; i < BULLETS_MAX; ++i) {
		if (!alive[i]) continue;
		p[i] += v2(0,1) * 300.0f * delta_time();
		if (p[i].y > 240.0f + 10.0f) {
			alive[i] = false;
		}
	}
}

void BulletBarn::draw()
{
	for (int i = 0; i < BULLETS_MAX; ++i) {
		if (!alive[i]) continue;
		v2 e = v2(2.0f, 5.0);
		draw_box(aabb(p[i] - e, p[i] + e), color_white());
	}
}

inline float flicker(float interval, float lo, float hi)
{
	if (between_interval(interval)) return lo;
	else return hi;
}

polygon make_asteroid_poly(v2 p0)
{
	polygon poly;
	for (int i = 0; i < POLYGON_MAX_VERTS; ++i) {
		float size = 30;
		poly.verts[i] = rnd_next_range(v2(-size,-size), v2(size,size));
		poly.verts[i] += p0;
	}
	poly.count = convex_hull(poly.verts, POLYGON_MAX_VERTS);
	poly.compute_norms();
	return poly;
}

void AsteroidBarn::add(v2 p, v2 v)
{
	for (int i = 0; i < ASTEROIDS_MAX; ++i) {
		if (alive[i]) continue;
		alive[i] = true;
		angular_velocity[i] = rnd_next_range(g->rnd, 0.2f, 1.0f) * sign(rnd_next_range(g->rnd, -1,1));
		velocity[i] = v;
		poly[i] = make_asteroid_poly(p);
		center_of_mass[i] = calc_center_of_mass(poly[i]);
		slice_timeout[i] = 0;
		return;
	}
}

void AsteroidBarn::add(polygon p, v2 v, float a, float timeout)
{
	for (int i = 0; i < ASTEROIDS_MAX; ++i) {
		if (alive[i]) continue;
		alive[i] = true;
		angular_velocity[i] = a;
		velocity[i] = v;
		poly[i] = p;
		center_of_mass[i] = calc_center_of_mass(poly[i]);
		slice_timeout[i] = timeout;
		return;
	}
}

void AsteroidBarn::update()
{
	for (int i = 0; i < ASTEROIDS_MAX; ++i) {
		if (!alive[i]) continue;

		// Delete tiny asteroids.
		if (calc_area(poly[i]) < 100.0f) {
			alive[i] = false;
			continue;
		}

		// Integrate positions.
		v2 delta = velocity[i] * delta_time();
		for (int j = 0; j < poly[i].count; ++j) {
			poly[i].verts[j] += delta;
		}
		center_of_mass[i] += delta;

		// Rotate about center of mass.
		// Integrate orientation.
		float angular_delta = angular_velocity[i] * delta_time();
		rotation r = sincos(angular_delta);
		v2 c = center_of_mass[i];
		for (int j = 0; j < poly[i].count; ++j) {
			poly[i].verts[j] = mul(r, poly[i].verts[j] - c) + c;
		}
		poly[i].compute_norms();

		if (slice_timeout[i] > 0) {
			slice_timeout[i] = max(0.0f, slice_timeout[i] - delta_time());
		}
	}
}

void AsteroidBarn::slice(ray r)
{
	halfspace h = halfspace(skew(r.d), r.p);
	for (int i = 0; i < ASTEROIDS_MAX; ++i) {
		if (!alive[i]) continue;
		if (slice_timeout[i] > 0) continue;
		if (!raycast(r, poly[i])) continue;
		sutherland_hodgman_output sho = sutherland_hodgman(h, poly[i]);
		if (!sho.front.count | !sho.back.count) continue;
		sho.front.compute_norms();
		sho.back.compute_norms();
		v2 c = center_of_mass[i];
		v2 v = velocity[i];
		float a = angular_velocity[i] * 2;
		v2 v_front = v + v2(-50.0f, 0);
		v2 v_back = v + v2(50.0f, 0);
		alive[i] = false;
		add(sho.front, v_front, a < 0 ? -a : a, 1.0f);
		add(sho.back, v_back, a < 0 ? a : -a, 1.0f);
	}
}

void AsteroidBarn::draw()
{
	for (int i = 0; i < ASTEROIDS_MAX; ++i) {
		if (!alive[i]) continue;
		draw_polygon(poly[i], color_white());
	}
}

void player_movement_routine()
{
	rt_begin(g->player.rt_movement, delta_time())
	{
		if (g->player.fired_laser) {
			nav_goto("laser_knockback");
		}
		v2 dir = v2(0,0);
		if (key_down('A')) {
			dir += v2(-1,0);
		}
		if (key_down('W')) {
			dir += v2(0,1);
		}
		if (key_down('D')) {
			dir += v2(1,0);
		}
		if (key_down('S')) {
			dir += v2(0,-1);
		}
		g->player.old_dir = g->player.dir;
		g->player.dir = dir;

		g->player.turn_time += delta_time();
		if (g->player.turn_time > 0.05f) {
			g->player.turn_time = 0;
			int target;
			if (g->player.dir.x == 0) {
				target = g->player.facing_index < 7 ? 3 : 10;
			} else if (g->player.dir.x > 0) {
				target = g->player.facing_index < 7 ? 6 : 13;
			} else {
				target = g->player.facing_index < 7 ? 0 : 7;
			}
			int diff = target - g->player.facing_index;
			int target_sign = diff != 0 ? sign(diff) : 0;
			g->player.facing_index = min(13, max(0, g->player.facing_index + target_sign));
			g->player.sprite.frame = g->player.facing_index;
		}

		bool go_slow = g->player.charging_laser | g->player.shielding | g->player.firing_rockets;
		g->player.p += safe_norm(dir) * (go_slow ? 100.0f : 200.0f) * delta_time();
		nav_restart();
	}

	rt_label("laser_knockback")
	{
		v2 offset = v2(0,-1) * 10.0f;
		int n = 5;
		for (int i = 0; i < n; ++i) {
			float t = (float)i / (float)n;
			v2 lo = lerp(t, g->player.bounds.min, g->player.bounds.min + offset);
			v2 hi = lerp(t, g->player.bounds.max, g->player.bounds.max + offset);
			g->trails.add(aabb(lo, hi), 0.15f);
		}
		g->player.p += offset;
	}
	rt_wait(0.25f)
	{
		nav_restart();
	}
	rt_end();

	g->player.bounds = aabb(g->player.p, 10, 20);

	g->player.booster_time += delta_time();
	if (g->player.booster_time > 0.025f) {
		g->player.booster_time = 0;
		if (g->player.facing_index < 7) {
			g->player.facing_index += 7;
		} else {
			g->player.facing_index -= 7;
		}
		g->player.sprite.frame = g->player.facing_index;
	}
}

void player_weapons_routine()
{
	rt_begin(g->player.rt_weapons, delta_time());
	{
		if (key_down('N')) {
			nav_goto("laser");
		} else if (key_pressed(TK_SPACE)) {
			nav_goto("bullet");
		} else if (key_pressed('M')) {
			nav_goto("rocket");
		} else if (key_pressed('B')) {
			nav_goto("shield");
		} else {
			nav_restart();
		}
	}

	rt_label("laser")
	{
		g->player.charging_laser = true;
	}
	// Charging laser.
	rt_seconds(1.0f)
	{
		if (key_down('N')) {
			v2 p = top(g->player.bounds) + v2(0, 10.0f);
			draw_circle(circle(p, fade(rt.elapsed, 3.0f, 10.0f)), color_white());
		} else {
			g->player.charging_laser = false;
			nav_restart();
		}
	}
	// Fire laser.
	rt_once()
	{
		v2 p = top(g->player.bounds) + v2(0, 10.0f);
		ray r = ray(p + v2(0,8), v2(0,1), 500);
		g->asteroids.slice(r);
		g->player.laser_trail = r;
		g->player.fired_laser = true;
	}
	// Laser Fade.
	rt_seconds(0.75f)
	{
		g->player.fired_laser = false;
		g->player.charging_laser = false;
		ray r = g->player.laser_trail;
		int alpha = fade_alpha(rt.elapsed, 0.4f, 1.0f);
		v2 p = r.p - v2(0,3);
		draw_line(p, r.endpoint(), color_white(alpha));
		draw_line(p + v2(1,0), p + v2(0,20) + v2(1,0), color_white(alpha));
		draw_line(p - v2(1,0), p + v2(0,20) - v2(1,0), color_white(alpha));
		draw_line(p + v2(2,0), p + v2(0,8) + v2(2,0), color_white(alpha));
		draw_line(p - v2(2,0), p + v2(0,8) - v2(2,0), color_white(alpha));
		draw_line(p + v2( 3,-1), p + v2(0,3) + v2( 3,-1), color_white(alpha));
		draw_line(p + v2(-3,-1), p + v2(0,3) + v2(-3,-1), color_white(alpha));
		draw_circle_fill(circle(r.p - v2(0,8), 5), color_white(alpha));
	}
	rt_once()
	{
		nav_restart();
	}

	rt_label("bullet")
	{
		g->player.bullets.add();
	}
	rt_seconds(0.75f)
	{
 		if (!key_down(TK_SPACE)) {
			nav_restart();
		}
	}
	rt_seconds(0.75f)
	{
		g->player.charging_shot = 1;
 		if (!key_down(TK_SPACE)) {
			nav_goto("bullet_cooldown");
		}
	}
	rt_seconds(key_down(TK_SPACE))
	{
		g->player.charging_shot = 2;
		if (key_down(TK_SPACE)) nav_redo();
		else {
			// Do charge shot.
			nav_restart();
		}
	}

	rt_label("bullet_cooldown") { }
	rt_wait(0.05f)
	{
		nav_restart();
	}

	rt_label("rocket")
	{
		g->player.firing_rockets = true;
	}
	rt_wait(0.1f)
	{
		v2 target = v2(g->player.p.x, 200) + rnd_next_range(v2(-40,-40), v2(20,20));
		g->player.rockets.add(g->player.p, target, 1.5f);
		if (g->player.rockets.count < g->player.rockets.capacity) {
			nav_goto("rocket");
		}
		g->player.firing_rockets = false;
	}
	rt_wait(1.0f)
	{
		nav_restart();
	}

	rt_label("shield")
	{
		g->player.shielding = true;
		g->player.shield_time += delta_time();
		if (key_down('B') && g->player.shield_time < player_shield_max) {
			nav_redo();
		}
		g->player.shield_time = 0;
		g->player.shielding = false;
	}
	rt_wait(0.5f)
	{
		nav_restart();
	}
	rt_end();
}

void set_global_pointers(dll_state* state)
{
	g_malloc = state->malloc;
	g_free = state->free;
	g = (game*)state->context;
	screen = state->tigr;
}

extern "C" __declspec( dllexport ) void game_init(dll_state* state)
{
	g_malloc = state->malloc;
	g_free = state->free;
	g = (game*)ALLOC(sizeof(game));
	memset(g, 0, sizeof(game));
	state->context = g;
	set_global_pointers(state);

	g->rnd = rnd_seed(0);
	g->player.reset();
	load_all_sprites();
}

extern "C" __declspec( dllexport ) void game_hotload(dll_state* state)
{
	set_global_pointers(state);
	memset(&g->asteroids, 0, sizeof(g->asteroids));
	g->player.reset();
}

extern "C" __declspec( dllexport ) void game_loop(float dt)
{
	update_time(dt);

	tigrClear(screen, color_black());

	player_movement_routine();
	player_weapons_routine();

	// Update + draw rockets.
	g->player.rockets.update();
	g->player.rockets.draw();
	v2 rocket_hit;
	while (g->player.rockets.try_pop_hit(&rocket_hit)) {
		g->player.rockets.count--;
	}

	// Update + draw bullets.
	g->player.bullets.update();
	g->player.bullets.draw();

	// Draw player.
	g->player.sprite.p = g->player.p;
	g->player.sprite.draw();
	if (g->player.shielding) {
		float shield_t = map(1.0f - (g->player.shield_time / player_shield_max), 0.4f, 1.0f);
		int shield_alpha = (int)(0.5f + shield_t * 255.0f);
		draw_circle(circle(g->player.p, flicker(0.075f, 25, 26)), color_white(shield_alpha));
	}

	// Update + draw asteroids.
	if (on_interval(0.3f)) {
		g->asteroids.add(g->player.p + v2(rnd_next_range(g->rnd,-50,50),600), v2(0,-100));
	}
	g->asteroids.update();
	g->asteroids.draw();

	// Update + draw trails.
	g->trails.update();
	g->trails.draw();

	g->anims.update();
	g->anims.draw();

	// Megaman buster for normal shot.
	// 3 stage (no charge, charged, max charged). 3rd stage has faster bullet speed + minor
	// knockback (on player when shot & against enemy on hit) + brief stun-shock (enemy).

	// Shield toggle.

	// Boss 1: shield (rotates)
	// --> beat normally w/ timing, stun with charge shot and spam minis, no stun lock
	// 3 max shots + barrage worth of hp.
	// 
	// Boss 2: laser (asteroids in between)
	// --> beat w/ shield reflect
	// 
	// Boss 3: missiles (hides behind asteroids/cover)
	// --> beat w/ laser penetration
	// 
	// Boss 4 (final?): fast + side-dash
	// --> beat w/ missile curve & AoE
}
