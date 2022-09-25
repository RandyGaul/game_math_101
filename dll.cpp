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

	int hit_count;
	v2 hits[ROCKETS_MAX];

	bool add(v2 start, v2 end, float duration);
	void update(float dt);
	void draw();
	bool try_pop_hit(v2* hit_out);
};

struct player_ship
{
	v2 p;
	aabb bounds;
	Coroutine co;
	rocket_barn rockets;

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
			p[i] = start;
			v2 d = norm(end - start);
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
	rockets.capacity = 3;
}

void rocket_barn::update(float dt)
{
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
		float t = elapsed[i] / hang_time[i];
		v2 at = bezier(p[i], c0[i], c1[i], target[i], ease_in_sin(t));
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
	g->player.p += safe_norm(dir) * 100.0f * dt;
	g->player.bounds = aabb(g->player.p, 10, 20);

	draw_box(g->player.bounds, color_white());

	co_begin(g->player.co, dt);
	co_for(1.0f)
	{
		if (tigrKeyHeld(screen, 'N')) {
			v2 p = top(g->player.bounds) + v2(0, 10.0f);
			draw_circle(circle(p, map(1.0f - smoothstep(co_for_t()), 3.0f, 10.0f)), color_white());
		} else {
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
	co_wait(1.0f);
	co_restart();
	co_end();

	if (tigrKeyDown(screen, 'M') && !tigrKeyHeld(screen, 'N')) {
		g->player.rockets.add(g->player.p, v2(g->player.p.x, 200), 1.5f);
	}

	g->player.rockets.update(dt);
	g->player.rockets.draw();
	v2 rocket_hit;
	while (g->player.rockets.try_pop_hit(&rocket_hit)) {
		g->player.rockets.count--;
	}
}
