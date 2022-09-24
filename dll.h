#pragma once

typedef void* (malloc_fn)(size_t);
typedef void (free_fn)(void*);

struct dll_state
{
	malloc_fn* malloc;
	free_fn* free;
	Tigr* tigr;
	void* context;
};

typedef void (game_init_fn)(dll_state*);
typedef void (game_hotload_fn)(dll_state*);
typedef void (game_loop_fn)(float dt);
