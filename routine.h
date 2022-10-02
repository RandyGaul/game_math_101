#ifndef ROUTINE_H
#define ROUTINE_H

#include <stdint.h>

// A portable "coroutine"-like thing for implementing FSM and behavior-cycles.
//
// Routine is a set of macros useful to implement Finite State Machine type of
// behaviors. It works really well for simpler behaviors with loops or cycles,
// and works especially well for anything resembling a cutscene.
// 
// The macros here are wrappers around a switch statement. All of the `rt_***`
// macros save a bookmark. On the next frame the routine is resumed at the last
// bookmark.
// 
// Each `rt_***` macro has a block { }. It contains code to run, and is attached
// to the `rt_***` macro. The exception is `rt_end` has no block. Blocks can have
// local variables but they don't persist, so be careful with them. A block will
// run for one frame.

struct Routine
{
	// A "hidden feature" - Goes from 0 to 1 over X seconds during `rt_seconds`.
	// Useful for animating things.
	float elapsed = 0;

	// Used in `rt_seconds` to repeat the block for X seconds.
	float seconds = 0;

	// Used in `rt_wait` to wait for X seconds.
    float wait_elapsed = 0;

	// Current block in the routine we're at.
	uint64_t at = 0;
};

// -------------------------------------------------------------------------------------------------
// Example:
// rt_begin(routine, delta_time())
// {
//     // This is the starting block. Runs once by default.
// }
// 
// rt_label("state 1")
// {
//     // Runs once.
// }
// rt_seconds(3)
// {
//     // Runs once per frame for 3 seconds.
// }
// rt_wait(1)
// {
//     // Waits 1 second, then runs this once.
// }
// rt_once()
// {
//     // Runs this one time.
//     
//     if (key_pressed(SPACE))
//         nav_goto("state 2"); // Next frame will begin on `rt_label("state 2")`.
//     else
//         nav_goto("state 1"); // Restart this state.
// }
// 
// rt_label("state 2") { } // Empty block.
// rt_wait(5.0f) { } // Empty block. This is like a long "cooldown" of 5 seconds.
// rt_once()
// {
//     nav_restart(); // Restart the whole routine back at `rt_begin`.0
// }
// rt_end();

// Begins the routine.
#define rt_begin(routine, dt)                                                \
    do {                                                                     \
        Routine& __rt = routine;                                             \
        bool __mn = true;                                                    \
        float __dt = dt;                                                     \
        if (__rt.wait_elapsed > 0) __rt.wait_elapsed -= dt;                  \
        else switch (__rt.at) {                                              \
        case 0: do {                                                         \

// Has a name that can be jumped to with `co_goto(name)`.
#define rt_label(name)                                                       \
        } while (0); if (__mn) __rt.at = rt_fnv1a(name); break;              \
        case rt_fnv1a(name): do {                                            \

// Runs its block once.
#define rt_once()                                                            \
        } while (0); if (__mn) __rt.at = __LINE__; break;                    \
        case __LINE__: do {                                                  \

// Runs its block for `time` seconds.
#define rt_seconds(time)                                                     \
        rt_once();                                                           \
        auto& rt = __rt;                                                     \
        if (__rt.seconds < time) {                                           \
            __rt.seconds += __dt;                                            \
            __mn = __rt.seconds >= time;                                     \
            if (__mn) {                                                      \
                __rt.elapsed = 1;                                            \
                __rt.seconds = 0;                                            \
            } else {                                                         \
                __rt.elapsed = __rt.seconds / (float)time;                   \
            }                                                                \
        }                                                                    \

// Runs its block while `condition` is true.
#define rt_while(condition)                                                  \
        } while (0); if (__mn) __rt.at = __LINE__; break;                    \
        case __LINE__: if (condition) {                                      \
            do { __mn = false;                                               \

// Runs its block once when `condition` becomes true.
#define rt_upon(condition)                                                   \
        } while (0); if (__mn) __rt.at = (condition) ? __LINE__ : -__LINE__; \
        break;                                                               \
        case -__LINE__:                                                      \
            if (condition) __rt.at = __LINE__;                               \
            break;                                                           \
        case __LINE__: do {                                                  \

// Waits for `time` seconds before running its block.
#define rt_wait(time)                                                        \
        } while (0); if (__mn) {                                             \
            __rt.wait_elapsed = time;                                        \
            __rt.at = __LINE__;                                              \
        }                                                                    \
        break;                                                               \
        case __LINE__: do {                                                  \

// End of the routine. Does not have a block.
#define rt_end()                                                             \
        } while (0); if (__mn) __rt.at = -1;                                 \
        break;                                                               \
        } __rt_end:;                                                         \
    } while (0);                                                             \

// -------------------------------------------------------------------------------------------------
// The `nav_***` macros can be used anywhere in a coroutine.
// Call them like normal functions (inside if-statements, loops, etc.).

// Repeats the block of code this is placed within.
// Skips the rest of the block.
#define nav_redo()                                                           \
        do { goto __rt_end; } while (0)                                      \

// Goes to the block with the name set by `rt_label(name)`.
#define nav_goto(name)                                                       \
        do {                                                                 \
            __rt.at = rt_fnv1a(name);                                        \
            goto __rt_end;                                                   \
        } while (0)                                                          \

// Restarts the whole routine.
#define nav_restart()                                                        \
        do {                                                                 \
            __rt.at = 0;                                                     \
            __rt.elapsed = 0;                                                \
            __rt.seconds = 0;                                                \
            __rt.wait_elapsed = 0;                                           \
            goto __rt_end;                                                   \
        } while (0)                                                          \

// Goes to `name`'s block on the next frame.
// Can be called *outside* of the coroutine.
#define nav_next(routine, name)                                              \
        do { routine.at = rt_fnv1a(name); } while (0)                        \

// -------------------------------------------------------------------------------------------------
// For internal use by `rt_label(name)` and `nav_goto(name)`.

inline uint64_t constexpr rt_fnv1a(const char* name)
{
	uint64_t h = 14695981039346656037ULL;
	char c = 0;
	while ((c = *name++)) {
		h = h ^ (uint64_t)c;
		h = h * 1099511628211ULL;
	}
	return h;
}

#endif // ROUTINE_H
