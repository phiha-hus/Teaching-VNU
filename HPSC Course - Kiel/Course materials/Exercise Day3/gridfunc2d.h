#ifndef GRIDFUNC2D_HEADER
#define GRIDFUNC2D_HEADER

#include "basic.h"
#include "linalg.h"
#include "miniblas.h"

struct _gridfunc2d
{
    uint n;   // Dimension of the gridfunction
    uint d;   // Dimension including the border. d == n + 2
    field *x; // Array of size d * d == (n+2) * (n+2).
};

typedef struct _gridfunc2d gridfunc2d;
typedef gridfunc2d *pgridfunc2d;

pgridfunc2d new_gridfunc2d(uint n);
void del_gridfunc2d(pgridfunc2d gf);

/**
 * Clears all values of a grid function 'gf', also the boundary values
 */
void clear_gridfunc2d(pgridfunc2d gf);

/**
 * Perform a vector addition.
 * 
 * y <-- y + alpha * x
 */
void add_gridfunc2d(pgridfunc2d y, field alpha, pgridfunc2d x);

/**
 * Write the contents of a grid function 'gf' to the file 'filename'
 */
void write_gridfunc2d(pgridfunc2d gf, const char *filename);

/**
 * Initialization of grid function with some sine function
 */
void init_sine_gridfunc2d(pgridfunc2d gf);

/**
 * Perform a matrix-vector multiplication of the 5-point-stencil L.
 * 
 * y <-- y + alpha * L * x
 */
void addeval_5point_stencil_gridfunc2d(pgridfunc2d y, field alpha, pgridfunc2d x);

/**
 * Black-box solver for the poisson problem with zero boundary condtions.
 * 
 * 'x' contains the right-hand-side on entry and contains the solution afterwards.
 * 
 * ATTENTION: the function 'addeval_5point_stencil_gridfunc2d' needs to be implemented
 * correctly in order to let the solver work correctly.
 */
void solve_gridfunc2d(pgridfunc2d x);

/**
 * Compute 2-norm of the difference of 'a' and 'b'.
 */
real diffnorm_gridfunc2d(pgridfunc2d a, pgridfunc2d b);

#endif //GRIDFUNC2D_HEADER