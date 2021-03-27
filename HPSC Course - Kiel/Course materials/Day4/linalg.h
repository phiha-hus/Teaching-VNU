#ifndef LINALG_HEADER
#define LINALG_HEADER

#include "basic.h"
#include "miniblas.h"

struct _vector
{
    field *x;
    uint dim;
};

typedef struct _vector vector;
typedef vector *pvector;

pvector new_vector(uint dim);
void del_vector(pvector x);
void clear_vector(pvector x);
void copy_vector(pvector src, pvector dest);
void random_vector(pvector x);
void add_vector(pvector x, field alpha, pvector y);

struct _matrix
{
    field *x;
    uint rows;
    uint cols;
    uint ld;
};

typedef struct _matrix matrix;
typedef matrix *pmatrix;

pmatrix new_matrix(uint rows, uint cols);
void del_matrix(pmatrix x);
void clear_matrix(pmatrix x);

/* ------------------------------------------------------------
 * UL decomposition using BLAS
 * ------------------------------------------------------------ */

void eval_l(const pmatrix a, pvector x);
void eval_u(const pmatrix a, pvector x);
void eval_ul(const pmatrix a, pvector x);
void decomp_ul(pmatrix a);

/* ------------------------------------------------------------
 * iterative solvers
 * ------------------------------------------------------------ */

/*
 * General type of a matrix-vector-multiplication.
 *
 * y <-- y + alpha A x,
 *
 * where A is given either implicitly or by 'matrix'
 */
typedef void (*addeval_func)(double alpha, void *matrix, pvector x, pvector y);

/*
 * CG-solver.
 *
 * 'A' and 'Matrix' define a matrix-vector-mulitplication. A has to be symmetric and positive-definite
 * 'x' is initially the right-hand-side of the linear system Ax=b to be solved.
 * After completion 'x' will contain the solution vector.
 */
void cg_solve(addeval_func A, void *matrix, pvector x);

#endif // LINALG_HEADER