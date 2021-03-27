#ifndef BASIC_HEADER
#define BASIC_HEADER

#include <stdio.h>
#include <stdlib.h>
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

#endif // BASIC_HEADER