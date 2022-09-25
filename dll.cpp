#include "tigr.h"
#include "math_101.h"
#include "draw.h"
#include "dll.h"
#include "coroutine.h"

#include <stdio.h>

malloc_fn* g_malloc;
free_fn* g_free;
#define ALLOC g_malloc
#define FREE g_free

#define ROCKETS_MAX 16

struct rocket_barn
{
	int count;
	int capacity;
	bool alive[ROCKETS_MAX];
	float elapsed[ROCKETS_MAX];
	float hang_time[ROCKETS_MAX];
	v2 p[ROCKETS_MAX];
	v2 c0[ROCKETS_MAX];
	v2 c1[ROCKETS_MAX];
	v2 target[ROCKETS_MAX];

	float prev_dt;
	int hit_count;
	v2 hits[ROCKETS_MAX];

	bool add(v2 start, v2 end, float duration);
	void update(float dt);
	void draw();
	bool try_pop_hit(v2* hit_out);
};

#define BULLETS_MAX 16

struct bullet_barn
{
	bool alive[BULLETS_MAX];
	v2 p[BULLETS_MAX];

	bool add();
	void update(float dt);
	void draw();
};

struct player_ship
{
	v2 p;
	aabb bounds;
	Coroutine co;
	bool firing_laser;
	bool firing_rockets;
	bool shielding;
	bool firing_bullet;
	float shield_time;
	rocket_barn rockets;
	bullet_barn bullets;

	void reset();
};

struct game
{
	rnd_t rnd;
	player_ship player;
};

game* g;

v2 rnd_next_range(v2 lo, v2 hi) { return v2(rnd_next_range(g->rnd, lo.x, hi.x), rnd_next_range(g->rnd, lo.y, hi.y)); }

bool rocket_barn::add(v2 start, v2 end, float duration)
{
	if (count >= capacity) return false;
	for (int i = 0; i < ROCKETS_MAX; ++i) {
		if (!alive[i]) {
			alive[i] = true;
			elapsed[i] = 0;
			hang_time[i] = duration;
			v2 d = norm(end - start);
			p[i] = start;
			c0[i] = (start - d * 100.0f + skew(d) * rnd_next_range(g->rnd, 25.0f, 150.0f) * sign(rnd_next_range(g->rnd, -1.0f, 1.0f)));
			c1[i] = rnd_next_range((start + end) * 0.5f - v2(50,50), (start + end) * 0.5f+ v2(50,50));
			target[i] = end;
			++count;
			return true;
		}
	}
	return false;
}

void player_ship::reset()
{
	memset(this, 0, sizeof(*this));
	rockets.capacity = 5;
}

void rocket_barn::update(float dt)
{
	prev_dt = dt;
	hit_count = 0;
	for (int i = 0; i < ROCKETS_MAX; ++i) {
		if (!alive[i]) continue;
		elapsed[i] += dt;
		if (elapsed[i] > hang_time[i]) {
			alive[i] = false;
			hits[hit_count++] = target[i];
		}
	}
}

void rocket_barn::draw()
{
	for (int i = 0; i < ROCKETS_MAX; ++i) {
		if (!alive[i]) continue;
		float t0 = (elapsed[i] - prev_dt * 2) / hang_time[i];
		float t1 = (elapsed[i] - prev_dt) / hang_time[i];
		float t = elapsed[i] / hang_time[i];
		v2 at0 = bezier(p[i], c0[i], c1[i], target[i], ease_in_sin(t0));
		v2 at1 = bezier(p[i], c0[i], c1[i], target[i], ease_in_sin(t1));
		v2 at = bezier(p[i], c0[i], c1[i], target[i], ease_in_sin(t));
		draw_circle_fill(circle(at0, 4), color_white(0x88));
		draw_circle_fill(circle(at1, 5), color_white(0xBB));
		draw_circle(circle(at, 3), color_white());
	}
}

bool rocket_barn::try_pop_hit(v2* hit_out)
{
	if (hit_count) {
		*hit_out = hits[--hit_count];
		return true;
	} else {
		return false;
	}
}

bool bullet_barn::add()
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

void bullet_barn::update(float dt)
{
	for (int i = 0; i < BULLETS_MAX; ++i) {
		if (!alive[i]) continue;
		p[i] += v2(0,1) * 300.0f * dt;
		if (p[i].y > 240.0f + 10.0f) {
			alive[i] = false;
		}
	}
}

void bullet_barn::draw()
{
	for (int i = 0; i < BULLETS_MAX; ++i) {
		if (!alive[i]) continue;
		v2 e = v2(2.0f, 5.0);
		draw_box(aabb(p[i] - e, p[i] + e), color_white());
	}
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
}

extern "C" __declspec( dllexport ) void game_hotload(dll_state* state)
{
	set_global_pointers(state);
	g->player.reset();
}

inline float flicker(float* t_ptr, float dt, float interval, float lo, float hi)
{
	float t = *t_ptr += dt;
	if (t < interval) return lo;
	else {
		if (t > interval * 2.0f) *t_ptr = 0;
		return hi;
	}
}

extern "C" __declspec( dllexport ) void game_loop(float dt)
{
	tigrClear(screen, color_black());

	// Player movement.
	v2 dir = v2(0,0);
	if (tigrKeyHeld(screen, 'A')) {
		dir += v2(-1,0);
	}
	if (tigrKeyHeld(screen, 'W')) {
		dir += v2(0,1);
	}
	if (tigrKeyHeld(screen, 'D')) {
		dir += v2(1,0);
	}
	if (tigrKeyHeld(screen, 'S')) {
		dir += v2(0,-1);
	}
	bool go_slow = g->player.firing_laser | g->player.shielding | g->player.firing_rockets;
	g->player.p += safe_norm(dir) * (go_slow ? 100.0f : 200.0f) * dt;
	g->player.bounds = aabb(g->player.p, 10, 20);

	// Player controls.
	const float shield_max = 3.0f;
	co_begin(g->player.co, dt);
	{
		if (tigrKeyHeld(screen, 'N')) {
			co_goto("laser");
		} else if (tigrKeyDown(screen, TK_SPACE)) {
			co_goto("bullet");
		} else if (tigrKeyDown(screen, 'M')) {
			co_goto("rocket");
		} else if (tigrKeyDown(screen, 'B')) {
			co_goto("shield");
		} else {
			co_restart();
		}
	}
	co_label("laser")
	{
		g->player.firing_laser = true;
		co_for(1.0f)
		{
			if (tigrKeyHeld(screen, 'N')) {
				v2 p = top(g->player.bounds) + v2(0, 10.0f);
				draw_circle(circle(p, map(1.0f - smoothstep(co_for_t()), 3.0f, 10.0f)), color_white());
			} else {
				g->player.firing_laser = false;
				co_restart();
			}
		}
		co_while(tigrKeyHeld(screen, 'N'))
		{
			v2 p = top(g->player.bounds) + v2(0, 10.0f);
			ray r = ray(p + v2(0,8), v2(0,1), 500);
			draw_line(r.p, r.endpoint(), color_white());
			static float t = 0;
			draw_circle_fill(circle(p, flicker(&t, dt, 0.05f, 5.0f, 6.0f)), color_white());
		}
		co_wait(1.0f) // 1 second cooldown.
		{
			g->player.firing_laser = false;
		}
		co_restart();
	}
	co_label("bullet")
	{
		g->player.firing_bullet = true;
		g->player.bullets.add();
		co_wait(0.05f);
		g->player.firing_bullet = false;
		co_restart();
	}
	co_label("rocket")
	{
		g->player.firing_rockets = true;
		g->player.rockets.add(g->player.p, v2(g->player.p.x, 200), 1.5f);
		co_wait(0.1f);
		if (g->player.rockets.count < g->player.rockets.capacity) {
			co_goto("rocket");
		}
		g->player.firing_rockets = false;
		co_wait(1.0f);
		co_restart();
	}
	co_label("shield")
	{
		g->player.shielding = true;
		g->player.shield_time += dt;
		if (tigrKeyHeld(screen, 'B') && g->player.shield_time < shield_max) {
			co_repeat();
		}
		g->player.shield_time = 0;
		g->player.shielding = false;
		co_wait(0.5f);
		co_restart();
	}
	co_end();

	// Update + draw rockets.
	g->player.rockets.update(dt);
	g->player.rockets.draw();
	v2 rocket_hit;
	while (g->player.rockets.try_pop_hit(&rocket_hit)) {
		g->player.rockets.count--;
	}

	// Update + draw bullets.
	g->player.bullets.update(dt);
	g->player.bullets.draw();

	// Draw player.
	draw_box(g->player.bounds, color_white());
	if (g->player.shielding) {
		float shield_t = map(1.0f - (g->player.shield_time / shield_max), 0.4f, 1.0f);
		int shield_alpha = (int)(0.5f + shield_t * 255.0f);
		static float t = 0;
		draw_circle(circle(g->player.p, flicker(&t, dt, 0.075f, 25, 26)), color_white(shield_alpha));
	}
}
