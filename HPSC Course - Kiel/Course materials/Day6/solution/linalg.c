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

pmatrix clone_matrix(pmatrix x)
{
    uint rows = x->rows;
    uint cols = x->cols;
    uint ld = x->ld;

    uint i, j;

    pmatrix y;

    y = new_matrix(rows, cols);

    for (j = 0; j < cols; j++)
    {
        for (i = 0; i < rows; i++)
        {
            y->x[i + j * ld] = x->x[i + j * ld];
        }
    }

    return y;
}

pbandmatrix new_bandmatrix(uint rows, uint cols, uint kl, uint ku)
{
    pbandmatrix x;

    x = (pbandmatrix)allocmem(sizeof(bandmatrix));
    x->rows = rows;
    x->cols = cols;
    x->kl = kl;
    x->ku = ku;
    x->x = (field *)allocmem((kl + ku + 1) * cols * sizeof(field));

    return x;
}
void del_bandmatrix(pbandmatrix x)
{
    assert(x != NULL);

    freemem(x->x);
    freemem(x);
}

void clear_bandmatrix(pbandmatrix x)
{
    uint rows = x->kl + x->ku + 1;
    uint cols = x->cols;
    field *xx = x->x;

    uint i, j;

    for (j = 0; j < cols; j++)
    {

        for (i = 0; i < rows; i++)
        {
            xx[i + j * rows] = 0.0;
        }
    }
}

pbandmatrix setup_5point_stencil_bandmatrix(uint n)
{
    uint ku = n;
    uint kl = n;
    uint rows = kl + ku + 1;
    uint cols = n * n;
    field *xx;

    uint i, j, off;

    pbandmatrix x;

    x = new_bandmatrix(n * n, n * n, n, n);
    xx = x->x;

    clear_bandmatrix(x);

    // put outermost superdiagonal -1
    for (j = n; j < cols; j++)
    {
        xx[j * rows] = -1.0;
    }

    // put first superdiagonal -1
    off = 1;
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n - 1; i++)
        {
            xx[(n - 1) + (off + i) * rows] = -1.0;
        }
        off += n;
    }

    // put main diagonal 4
    for (j = 0; j < cols; j++)
    {
        xx[n + j * rows] = 4.0;
    }

    // put first subdiagonal -1
    off = 0;
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n - 1; i++)
        {
            xx[(n + 1) + (off + i) * rows] = -1.0;
        }
        off += n;
    }

    // put outermost subdiagonal -1
    for (j = 0; j < cols - n; j++)
    {
        xx[(rows - 1) + j * rows] = -1.0;
    }

    return x;
}

void print_bandstorage_bandmatrix(pbandmatrix x)
{
    uint kl = x->kl;
    uint ku = x->ku;
    uint rows = kl + ku + 1;
    uint cols = x->cols;
    field *xx = x->x;

    uint i, j;

    assert(x != NULL);

    printf("Bandmatrix: %d rows, %d cols, %d bands:\n", x->rows, x->cols, rows);
    for (i = 0; i < rows; i++)
    {
        printf("(");
        if (i < ku)
        {
            for (j = 0; j < ku - i; j++)
            {
                printf("    *     ");
            }
            for (j = ku - i; j < cols; j++)
            {
                printf("%+.2e ", xx[i + j * rows]);
            }
        }
        else if (i == ku)
        {
            for (j = 0; j < cols; j++)
            {
                printf("%+.2e ", xx[i + j * rows]);
            }
        }
        else
        {
            for (j = 0; j < cols - (i - ku); j++)
            {
                printf("%+.2e ", xx[i + j * rows]);
            }

            for (j = cols - (i - ku); j < cols; j++)
            {
                printf("    *     ");
            }
        }
        printf(")\n");
    }
}

void decomplr_bandmatrix(pbandmatrix x)
{
    int kl = x->kl;
    int ku = x->ku;
    int bands = kl + ku + 1;
    int n = x->cols;
    field *xx = x->x;

    int i;
    for (i = 0; i + 1 < n; i++)
    {
        scal(UINT_MIN(kl, n - i - 1), 1.0 / xx[ku + i * bands], xx + (ku + 1) + i * bands, 1);

        ger(UINT_MIN(kl, n - i - 1), UINT_MIN(ku, n - i - 1), -1.0,
            xx + (ku + 1) + i * bands, 1,
            xx + (ku - 1) + (i + 1) * bands, bands - 1,
            xx + ku + (i + 1) * bands, bands - 1);
    }
}

void leval(int n, field *a, int lda, field *x, int incx)
{
    int k;

    for (k = n; k-- > 0;)
    {
        x[k] += dot(k, a + k, lda, x, incx);
    }
}

void reval(int n, field *a, int lda, field *x, int incx)
{
    int k;

    for (k = 0; k < n; k++)
    {
        x[k] = dot(n - k, a + k + (k)*lda, lda, x + k, incx);
    }
}

void rleval(int n, field *a, int lda, field *x, int incx)
{
    leval(n, a, lda, x, incx);
    reval(n, a, lda, x, incx);
}

void rldecomp(int n, field *a, int lda)
{
    int k;

    for (k = n; k-- > 0;)
    {
        scal(k, 1.0 / a[k + k * lda], a + k, lda);
        ger(k, k, -1.0, a + k * lda, 1, a + k, lda, a, lda);
    }
}

void block_rsolve(int n, int m, const field *r, int ldr,
                  field *x, int ldx)
{
    int k;

    for (k = n; k-- > 0;)
    {
        scal(m, 1.0 / r[k + k * ldr], x + k, ldx);
        ger(k, m, -1.0, r + k * ldr, 1, x + k, ldx, x, ldx);
    }
}

void block_lsolve_trans(int n, int m, const field *l, int ldl,
                        field *x, int ldx)
{
    int k;

    for (k = n; k-- > 0;)
    {
        ger(m, k, -1.0, x + k * ldx, 1, l + k, ldl, x, ldx);
    }
}

void block_rldecomp(int n, field *a, int lda, uint blocks)
{
    int i, j, k;
    int oi, oj, ok, ni, nj, nk;

    for (k = blocks; k-- > 0;)
    {
        ok = n * k / blocks;
        nk = n * (k + 1) / blocks - ok;
        rldecomp(nk, a + ok + ok * lda, lda);

        for (j = 0; j < k; j++)
        {
            oj = n * j / blocks;
            nj = n * (j + 1) / blocks - oj;
            block_rsolve(nk, nj,
                         a + ok + ok * lda, lda,
                         a + ok + oj * lda, lda);
            block_lsolve_trans(nk, nj,
                               a + ok + ok * lda, lda,
                               a + oj + ok * lda, lda);
        }

        for (j = 0; j < k; j++)
        {
            oj = n * j / blocks;
            nj = n * (j + 1) / blocks - oj;
            for (i = 0; i < k; i++)
            {
                oi = n * i / blocks;
                ni = n * (i + 1) / blocks - oi;

                gemm(false, false, ni, nj, nk, -1.0,
                     a + oi + ok * lda, lda,
                     a + ok + oj * lda, lda,
                     1.0, a + oi + oj * lda, lda);
            }
        }
    }
}

void par_block_rldecomp(int n, field *a, int lda, int blocks)
{
    int i, j, k;
    int oi, oj, ok, ni, nj, nk;

#pragma omp parallel
#pragma omp single
    {

        for (k = blocks; k-- > 0;)
        {
            ok = n * k / blocks;
            nk = n * (k + 1) / blocks - ok;
#pragma omp task depend(inout                \
                        : a[ok + ok * lda]), \
    firstprivate(nk, ok)
            rldecomp(nk, a + ok + ok * lda, lda);

            for (j = 0; j < k; j++)
            {
                oj = n * j / blocks;
                nj = n * (j + 1) / blocks - oj;
#pragma omp task depend(in                   \
                        : a[ok + ok * lda]), \
    depend(inout                             \
           : a[ok + oj * lda]),              \
    firstprivate(nk, nj, ok, oj)
                block_rsolve(nk, nj,
                             a + ok + ok * lda, lda,
                             a + ok + oj * lda, lda);
#pragma omp task depend(in                   \
                        : a[ok + ok * lda]), \
    depend(inout                             \
           : a[oj + ok * lda]),              \
    firstprivate(nk, nj, ok, oj)
                block_lsolve_trans(nk, nj,
                                   a + ok + ok * lda, lda,
                                   a + oj + ok * lda, lda);
            }

            for (j = 0; j < k; j++)
            {
                oj = n * j / blocks;
                nj = n * (j + 1) / blocks - oj;
                for (i = 0; i < k; i++)
                {
                    oi = n * i / blocks;
                    ni = n * (i + 1) / blocks - oi;

#pragma omp task depend(in                   \
                        : a[oi + ok * lda],  \
                          a[ok + oj * lda]), \
    depend(inout                             \
           : a[oi + oj * lda]),              \
    firstprivate(nk, nj, ni, ok, oj, oi)
                    gemm(false, false, ni, nj, nk, -1.0,
                         a + oi + ok * lda, lda,
                         a + ok + oj * lda, lda,
                         1.0, a + oi + oj * lda, lda);
                }
            }
        }
    }
}
