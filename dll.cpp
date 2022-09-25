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

struct player_ship
{
	v2 p;
	aabb bounds;
	Coroutine co;
};

struct game
{
	rnd_t rnd;
	player_ship player;
};

game* g;

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
}

extern "C" __declspec( dllexport ) void game_hotload(dll_state* state)
{
	set_global_pointers(state);
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
		if (tigrKeyHeld(screen, TK_SPACE)) {
			v2 p = top(g->player.bounds) + v2(0, 10.0f);
			draw_circle(circle(p, map(1.0f - smoothstep(co_for_t()), 3.0f, 10.0f)), color_white());
		} else {
			co_restart();
		}
	}
	co_while(tigrKeyHeld(screen, TK_SPACE))
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
}
