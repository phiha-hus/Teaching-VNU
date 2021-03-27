#include "linalg.h"

pvector new_vector(uint dim)
{
    pvector x;

    x = allocmem(sizeof(vector));
    x->dim = dim;
    x->x = allocmem(dim * sizeof(field));

    return x;
}

void del_vector(pvector x)
{
    assert(x != NULL);

    freemem(x->x);
    freemem(x);
}

void clear_vector(pvector x)
{
    uint dim = x->dim;

    uint i;
    for (i = 0; i < dim; i++)
    {
        x->x[i] = 0.0;
    }
}

void copy_vector(pvector src, pvector dest)
{
    int n = src->dim;

    int i;

    assert(dest->dim == n);

    for (i = 0; i < n; i++)
    {
        dest->x[i] = src->x[i];
    }
}

void random_vector(pvector x)
{
    double *xx = x->x;
    int rows = x->dim;
    int i;

    for (i = 0; i < rows; i++)
    {
        xx[i] = 2.0 * rand() / RAND_MAX - 1.0;
    }
}

void add_vector(pvector x, field alpha, pvector y)
{
    uint dim = x->dim;

    uint i;

    assert(dim == y->dim);

    for (i = 0; i < dim; i++)
    {
        x->x[i] += alpha * y->x[i];
    }
}

pmatrix new_matrix(uint rows, uint cols)
{
    pmatrix x;

    x = allocmem(sizeof(matrix));
    x->rows = rows;
    x->cols = cols;
    x->ld = rows;
    x->x = allocmem(rows * cols * sizeof(field));

    return x;
}

void del_matrix(pmatrix x)
{

    assert(x != NULL);

    freemem(x->x);
    freemem(x);
}

void clear_matrix(pmatrix x)
{
    uint rows = x->rows;
    uint cols = x->cols;
    uint ld = x->ld;

    uint i, j;

    for (j = 0; j < cols; j++)
    {
        for (i = 0; i < rows; i++)
        {
            x->x[i + j * ld] = 0.0;
        }
    }
}

/* ------------------------------------------------------------
 * UL decomposition using BLAS
 * ------------------------------------------------------------ */

void eval_l(const pmatrix a, pvector x)
{
    int n = a->rows;
    int lda = a->ld;
    double *aa = a->x;
    double *xx = x->x;
    int k;

    assert(a->cols == n);
    assert(x->dim == n);

    for (k = n; k-- > 0;)
    {
        axpy(n - k - 1, xx[k], aa + k + 1 + k * lda, 1, xx + k + 1, 1);
        //     xx[k] += dot(k, aa+k, lda, xx, 1);
    }
}

void eval_u(const pmatrix a, pvector x)
{
    int n = a->rows;
    int lda = a->ld;
    double *aa = a->x;
    double *xx = x->x;
    int k;

    assert(a->cols == n);
    assert(x->dim == n);

    for (k = 0; k < n; k++)
    {
        axpy(k, xx[k], aa + k * lda, 1, xx, 1);
        xx[k] *= aa[k + k * lda];
        //     xx[k] = dot(n-k, aa+k+(k)*lda, lda, xx+k, 1);
    }
}

void eval_ul(const pmatrix a, pvector x)
{
    eval_l(a, x);
    eval_u(a, x);
}

void decomp_ul(pmatrix a)
{
    int n = a->rows;
    int lda = a->ld;
    double *aa = a->x;
    int k;

    assert(a->cols == n);

    for (k = n; k-- > 0;)
    {
        scal(k, 1.0 / aa[k + k * lda], aa + k, lda);
        ger(k, k, -1.0,
            aa + k * lda, 1,
            aa + k, lda,
            aa, lda);
    }
}

#define EPS_SOLVE 1.0e-12
#define MAX_STEPS 1000000

void cg_solve(addeval_func A, void *matrix, pvector x)
{
    /**
     * TODO:
     */
}