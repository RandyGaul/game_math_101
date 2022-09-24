
#include "tigr.h"
#include "math_101.h"
#include "draw.h"
#include "dll.h"

malloc_fn* g_malloc;
free_fn* g_free;
#define ALLOC g_malloc
#define FREE g_free

struct game
{
	rnd_t rnd;
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
	state->context = g;
	set_global_pointers(state);

	g->rnd = rnd_seed(0);
}

extern "C" __declspec( dllexport ) void game_hotload(dll_state* state)
{
	set_global_pointers(state);
}

extern "C" __declspec( dllexport ) void game_loop(float dt)
{
	tigrClear(screen, color_red());
}
