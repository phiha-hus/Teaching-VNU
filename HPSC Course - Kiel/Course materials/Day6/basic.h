#ifndef BASIC_HEADER
#define BASIC_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include <omp.h>

/**
 * Define basic datatypes like real number, intergers.
 */

typedef double real;
typedef double field;

typedef unsigned int uint;

/** @brief Stop watch for run-time measurements. */
typedef struct _stopwatch stopwatch;

/** @brief Pointer to a @ref stopwatch object. */
typedef stopwatch *pstopwatch;

/** @brief Create a @ref stopwatch object.
 *
 *  @returns New @ref stopwatch object. */
pstopwatch
new_stopwatch();

/** @brief Delete a @ref stopwatch object.
 *
 *  @param sw Stopwatch object to be deleted. */
void del_stopwatch(pstopwatch sw);

/** @brief Start a stopwatch.
 *
 *  This function stores the current time in a private
 *  variable that can be used subsequently to measure elapsed time.
 *
 *  @param sw Stopwatch to be started. */
void start_stopwatch(pstopwatch sw);

/** @brief Stop a stopwatch.
 *
 *  This function stores the current time in a private
 *  variable that can be used subsequently to determine the time
 *  that has passed between calls to @ref start_stopwatch and
 *  @ref stop_stopwatch.
 *
 *  @param sw Stopwatch (@ref start_stopwatch has to have been called for this object at least once).
 *  @returns Elapsed time in seconds since stopwatch was started. */
real stop_stopwatch(pstopwatch sw);

/**
 * Memory allocation / deallocation
 */

#define allocmem(size) malloc(size)
#define freemem(pointer) free(pointer)

#define UINT_MIN(x, y) ((x) < (y) ? (x) : (y))
#define UINT_MAX(x, y) ((x) < (y) ? (y) : (x))

#define ABS(x) fabs(x)
#define SQRT(x) sqrt(x)

real sqr(real x);

#endif // BASIC_HEADER
