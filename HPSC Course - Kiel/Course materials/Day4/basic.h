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

/**
 * Memory allocation / deallocation
 */

#define allocmem(size) malloc(size)
#define freemem(pointer) free(pointer)

/**
 * Little helpers
 */

#define ABS(x) (((x) < 0.0) ? (-(x)) : (x))

#endif // BASIC_HEADER