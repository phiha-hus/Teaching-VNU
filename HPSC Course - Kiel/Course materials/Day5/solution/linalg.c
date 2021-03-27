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
    int i, n;
    pvector r, p, a;
    double norm, error, gamma, lambda, mu;

    n = x->dim;
    r = new_vector(n);
    p = new_vector(n);
    a = new_vector(n);

    norm = nrm2(n, x->x, 1); // Norm of rhs

    // Init residual vector 'r'
    scal(n, 0.0, r->x, 1);          // set r <-- 0
    axpy(n, 1.0, x->x, 1, r->x, 1); // copy x --> r
    A(-1.0, matrix, x, r);          // r <-- b-Ax

    // Init search direction 'p'
    scal(n, 0.0, p->x, 1);          // set p <-- 0
    axpy(n, 1.0, r->x, 1, p->x, 1); // copy r --> p

    // Compute error
    error = nrm2(n, r->x, 1);

    i = 0;
    while (error > EPS_SOLVE * norm && i + 1 != MAX_STEPS)
    {
        scal(n, 0.0, a->x, 1); // set a <-- 0
        A(1.0, matrix, p, a);  // a = A p

        gamma = dot(n, p->x, 1, a->x, 1); // lambda = <p, r> / <p, a>
        lambda = dot(n, p->x, 1, r->x, 1) / gamma;

        axpy(n, lambda, p->x, 1, x->x, 1);  // x = x + lambda p
        axpy(n, -lambda, a->x, 1, r->x, 1); // r = r - lambda a

        mu = dot(n, r->x, 1, a->x, 1) / gamma; // p = r - mu p
        scal(n, -mu, p->x, 1);
        axpy(n, 1.0, r->x, 1, p->x, 1);

        // Compute error
        error = nrm2(n, r->x, 1);

        i++;
    }

    del_vector(r);
    del_vector(p);
    del_vector(a);
}