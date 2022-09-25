#ifndef COROUTINE_H
#define COROUTINE_H

#include <stdint.h>

// Original implementation from Noel Berry.

struct Coroutine
{
	// Current time to wait until we run the next block.
	float wait_for = 0;

	// Used in `co_for` to repeat the block for X time.
	float repeat_for = 0;

	// Used in `co_for_t` to calculate a number from [0-1).
	float repeat_for_t = 0;

	// Current block in the coroutine we're at.
	uint64_t at = 0;
};

inline uint64_t constexpr co_fnv1a(const char* name)
{
	uint64_t h = (uint64_t)14695981039346656037U;
	char c = 0;
	while ((c = *name++)) {
		h = h ^ (uint64_t)c;
		h = h * (uint64_t)1099511628211;
	}
	return h;
}

// Begins the coroutine.
#define co_begin(coroutine, dt)                                              \
    do {                                                                     \
        Coroutine& __co = coroutine;                                         \
        bool __mn = true;                                                    \
        float __dt = dt;                                                     \
        if (__co.wait_for > 0) __co.wait_for -= dt;                          \
        else switch (__co.at) {                                              \
        case 0: {                                                            \

// Waits until the next frame and runs the following block once.
#define co_step()                                                            \
        } if (__mn) __co.at = __LINE__; break;                               \
        case __LINE__: {                                                     \

// Waits until the next frame to run the following block once.
// Has a name that can be jumped to with `co_goto(name)`.
#define co_label(name)                                                       \
        } if (__mn) __co.at = co_fnv1a(name); break;                         \
        case co_fnv1a(name): {                                               \

// Repeats the following block for the given time.
#define co_for(time)                                                         \
        co_step();                                                           \
        __co.repeat_for_t = min(__co.repeat_for / time, 1.0f);               \
        if (__co.repeat_for < time) {                                        \
            __co.repeat_for += __dt;                                         \
            __mn = __co.repeat_for >= time;                                  \
            if (__mn) __co.repeat_for = 0;                                   \
        }                                                                    \

// Returns a number from [0-1) during `co_for`.
// Returns 0 as `co_far` starts, and approaches 1 as `co_for` ends.
#define co_for_t()                                                           \
        (__co.repeat_for_t)                                                  \

// Repeats the following block while the condition is met.
#define co_while(condition)                                                  \
        } if (__mn) __co.at = __LINE__; break;                               \
        case __LINE__: if (condition) {                                      \
            __mn = false;                                                    \

// Waits until the condition is meet to run the following block.
#define co_until(condition)                                                  \
        } if (__mn) __co.at = (condition) ? __LINE__ : -__LINE__;            \
        break;                                                               \
        case -__LINE__:                                                      \
            if (condition) __co.at = __LINE__;                               \
            break;                                                           \
        case __LINE__: {                                                     \

// Waits a given amount of time before beginning the follow block.
#define co_wait(time)                                                        \
        } if (__mn) {                                                        \
            __co.wait_for = time;                                            \
            __co.at = __LINE__;                                              \
        }                                                                    \
        break;                                                               \
        case __LINE__: {                                                     \

#define co_end()                                                             \
        } if (__mn) __co.at = -1;                                            \
        break;                                                               \
        }                                                                    \
    } while (0); __co_end:                                                   \

// Repeats the block of code this is placed within.
// Skips the rest of the block.
#define co_repeat()                                                          \
        do { goto __co_end; } while (0)                                      \

// Goes to the block with the given name set by `co_label`.
#define co_goto(name)                                                        \
        do {                                                                 \
            __co.at = co_fnv1a(name);                                        \
            goto __co_end;                                                   \
        } while (0)                                                          \

// Restarts the whole coroutine.
// Jumps back to `co_begin` on the next frame.
#define co_restart()                                                         \
        do {                                                                 \
            __co.at = 0;                                                     \
            __co.wait_for = 0;                                               \
            __co.repeat_for = 0;                                             \
            __co.repeat_for_t = 0;                                           \
            goto __co_end;                                                   \
        } while (0)                                                          \

#endif // COROUTINE_H
