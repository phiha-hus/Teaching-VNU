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

typedef void (*addeval_t)(pgridfunc2d y, field alpha, pgridfunc2d x, void *data);

pgridfunc2d new_gridfunc2d(uint n);
void del_gridfunc2d(pgridfunc2d gf);

void clear_gridfunc2d(pgridfunc2d gf);

/**
 * Perform a vector addition.
 * 
 * y <-- y + alpha * x
 */
void add_gridfunc2d(pgridfunc2d y, field alpha, pgridfunc2d x);

/**
 * Performs the dot product of two grid function 'x' and 'y'
 */
field dot_gridfunc2d(pgridfunc2d x, pgridfunc2d y);

/**
 * Computes the 2-norm of a grid function 'x'
 */
real norm2_gridfunc2d(pgridfunc2d x);

/**
 * Scales a grid function 'x' by the scalar 'alpha'.
 * x <-- alpha *x
 */
void scale_gridfunc2d(pgridfunc2d x, field alpha);

void write_gridfunc2d(pgridfunc2d gf, const char *filename);

void init_sine_gridfunc2d(pgridfunc2d gf);

/**
 * Perform a matrix-vector multiplication of the 5-point-stencil L.
 * 
 * y <-- y + alpha * L * x
 */
void addeval_5point_stencil_gridfunc2d(pgridfunc2d y, field alpha, pgridfunc2d x, void *data);

/**
 * Evaluate the system matrix of the Crank-Nicolson scheme for the 2D wave equation
 */
void addeval_cn_wave2d(pgridfunc2d y, field alpha, pgridfunc2d x, void *data);

/**
 * Solve a linear system of equation by the conjugant gradient method
 * 
 * 'A' callback function for the matrix-vector product
 * 'x' Initially contains the right-hand-side of the linear system.
 *     After completion contains the solution of the system.
 * 'data' Additional data passed to 'A'.
 */
void solve_cg_gridfunc2d(addeval_t A, pgridfunc2d x, void *data);

#endif //GRIDFUNC2D_HEADER