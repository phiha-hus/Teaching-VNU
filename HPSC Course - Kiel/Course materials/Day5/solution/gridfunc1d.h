#ifndef GRIDFUNC1D_HEADER
#define GRIDFUNC1D_HEADER

#include "basic.h"
#include "linalg.h"
#include "miniblas.h"

struct _gridfunc1d
{
    uint n;   // Dimension of the gridfunction
    uint d;   // Dimension including the border d == n + 2
    field *x; // Array of size d == n + 2
};

typedef struct _gridfunc1d gridfunc1d;
typedef gridfunc1d *pgridfunc1d;

pgridfunc1d new_gridfunc1d(uint n);
void del_gridfunc1d(pgridfunc1d gf);

/**
 * Clears all values of a grid function 'gf', also the boundary values
 */
void clear_gridfunc1d(pgridfunc1d gf);

/**
 * Perform a vector addition.
 * 
 * y <-- y + alpha * x
 */
void add_gridfunc1d(pgridfunc1d y, field alpha, pgridfunc1d x);

/**
 * Write the contents of a grid function 'gf' to the file 'filename'
 */
void write_gridfunc1d(pgridfunc1d gf, const char *filename);

/**
 * Initialization of grid function with some sine function
 */
void init_sine_gridfunc1d(pgridfunc1d gf);

/**
 * Perform a matrix-vector multiplication of the 3-point-stencil.
 * 
 * y <-- y + alpha * L * x
 */
void addeval_3point_stencil_gridfunc1d(pgridfunc1d y, field alpha, pgridfunc1d x);

/**
 * CG solver for the poisson problem with zero boundary condtions.
 * 
 * 'x' contains the right-hand-side on entry and contains the solution afterwards.
 * 
 */
void cg_gridfunc1d(pgridfunc1d x);

/**
 * One step of the Richardson iteration for the poisson problem with zero boundary condtions.
 * 
 * 'theta' is a damping parameter and has to be set based on the problem dimension,
 * the larger the problem dimension, the smaller it has to be
 * 'x' contains the starting value on entry and contains the result afterwards.
 * 'b' contains the right-hand-side on entry and is not changed.
 * 'd' is an auxiliary vector.
 * 
 */
void rich_step_1d(real theta, pgridfunc1d b, pgridfunc1d x, pgridfunc1d d);

/**
 * One step of the Jacobi iteration for the poisson problem with zero boundary condtions.
 * 
 * 'theta' is a damping parameter and can in our case be set to around 0.5.
 * 'x' contains the starting value on entry and contains the result afterwards.
 * 'b' contains the right-hand-side on entry and is not changed.
 * 'd' is an auxiliary vector.
 * 
 */
void jac_step_1d(real theta, pgridfunc1d b, pgridfunc1d x, pgridfunc1d d);

/**
 * One step of the Gauß-Seidel iteration for the poisson problem with zero boundary condtions.
 * 
 * 'x' contains the starting value on entry and contains the result afterwards.
 * 'b' contains the right-hand-side on entry and is not changed.
 * 
 */
void gs_step_1d(pgridfunc1d b, pgridfunc1d x);

/**
 * Restricting a gridfunction from one grid to the next coarser grid.
 * 
 * 'f' is the gridfunction on the fine grid.
 * 'c' is the gridfunction on the coarse grid.
 * 
 *  c <-- r(f)
 */
void restrict_gridfunc1d(pgridfunc1d f, pgridfunc1d c);

/**
 * Prolonging a gridfunction from one grid to the next finer grid and adding it to a function on a the fine grid.
 * 
 * 'c' is the gridfunction on the coarse grid.
 * 'f' is the gridfunction on the fine grid.
 * 
 *  f <-- f + p(c)
 */
void prolong_gridfunc1d(pgridfunc1d c, pgridfunc1d f);


/**
 * Multigrid algorithm for the 1D model problem on a series of (L+1) nested grids with the
 * coarsest grid being 0 and the finest L
 * 
 * 'x' contains an array of gridfunctions corresponding to the grid hierarchy with x[L] being the 
 *  starting guess on the finest grid on entry and the solution on the finest grid afterwards.
 * 'b' contains an array of gridfunctions corresponding to the grid hierarchy with b[L] being the 
 *  right-hand side on the finest grid on entry.
 * 'd' contains an auxiliary array of gridfunctions corresponding to the grid hierarchy
 *
 */
void mg_gridfunc1d(uint L, pgridfunc1d *x, pgridfunc1d *b, pgridfunc1d *d);



/**
 * Compute 2-norm of the difference of 'a' and 'b'.
 */
real diffnorm_gridfunc1d(pgridfunc1d a, pgridfunc1d b);








#endif //GRIDFUNC1D_HEADER
