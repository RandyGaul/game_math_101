#define _CRT_NONSTDC_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <Windows.h>
#include <stdio.h>
#include <stdlib.h>

#include "tigr.h"
#include "dll.h"

dll_state state;
game_init_fn* game_init;
game_hotload_fn* game_hotload;
game_loop_fn* game_loop;
HMODULE game_dll;
FILETIME game_dll_last_write_time;

FILETIME get_last_write_time(const char* path)
{
	FILETIME time = { 0 };
	WIN32_FILE_ATTRIBUTE_DATA data;
	if (GetFileAttributesExA(path, GetFileExInfoStandard, &data)) time = data.ftLastWriteTime;
	return time;
}

void unload_dll()
{
	FreeLibrary(game_dll);
	game_dll = 0;
	game_loop = 0;
}

void reload_dll()
{
	unload_dll();

	WIN32_FILE_ATTRIBUTE_DATA unused;
	if (!GetFileAttributesExA("lock.tmp", GetFileExInfoStandard, &unused)) {
		CopyFileA("bin/game.dll", "bin/game_temp.dll", 0);
		game_dll = LoadLibraryA("bin/game_temp.dll");

		if (!game_dll) {
			DWORD err = GetLastError();
			printf("Can't load lib: %d\n", err);
			return;
		}

		game_init = (game_init_fn*)GetProcAddress(game_dll, "game_init");
		if (!game_init) {
			DWORD err = GetLastError();
			printf("Cant load game_init: %d\n", err);
			return;
		}

		game_hotload = (game_hotload_fn*)GetProcAddress(game_dll, "game_hotload");
		if (!game_hotload) {
			DWORD err = GetLastError();
			printf("Cant load game_hotload: %d\n", err);
			return;
		}

		game_loop = (game_loop_fn*)GetProcAddress(game_dll, "game_loop");
		if (!game_loop) {
			DWORD err = GetLastError();
			printf("Cant load game_loop: %d\n", err);
			return;
		}

		game_dll_last_write_time = get_last_write_time("bin/game_temp.dll");
	}
}

bool reload_dll_if_needed()
{
	FILETIME new_time = get_last_write_time("bin/game.dll");
	if (CompareFileTime(&new_time, &game_dll_last_write_time)) {
		reload_dll();
		return true;
	} else {
		return false;
	}
}

int main()
{
	Tigr* screen = tigrWindow(640, 480, "Math 101", 0);

	reload_dll();
	state.tigr = screen;
	state.malloc = malloc;
	state.free = free;
	game_init(&state);

	while (!tigrClosed(screen) && !tigrKeyDown(screen, TK_ESCAPE)) {
		float dt = tigrTime();
		if (tigrKeyHeld(screen, TK_CONTROL) && (tigrKeyDown(screen, TK_F5) || tigrKeyDown(screen, 'R'))) {
			system("build.cmd");
		}else if (reload_dll_if_needed()) {
			game_hotload(&state);
		}
		game_loop(dt);
		tigrUpdate(screen);
	}

	tigrFree(screen);

	return 0;
}
