#include "gridfunc2d.h"

pgridfunc2d new_gridfunc2d(uint n)
{
    pgridfunc2d gf;

    gf = (pgridfunc2d)allocmem(sizeof(gridfunc2d));
    gf->n = n;
    gf->d = n + 2;
    gf->x = (field *)allocmem(gf->d * gf->d * sizeof(field));

    return gf;
}

void del_gridfunc2d(pgridfunc2d gf)
{
    assert(gf != NULL);

    freemem(gf->x);
    freemem(gf);
}

void clear_gridfunc2d(pgridfunc2d gf)
{
    field *x;
    uint d;
    uint i, j;

    assert(gf != NULL);
    x = gf->x;
    d = gf->d;

    for (j = 0; j < d; j++)
    {
        for (i = 0; i < d; i++)
        {
            x[i + j * d] = 0.0;
        }
    }
}

void add_gridfunc2d(pgridfunc2d y, field alpha, pgridfunc2d x)
{
    field *xx, *yx;
    uint d, n;
    uint i, j;

    assert(y != NULL);
    assert(x != NULL);

    xx = x->x;
    yx = y->x;
    n = y->n;
    d = y->d;

    assert(x->d == d);
    assert(x->n == n);

    for (j = 1; j < n; j++)
    {
        for (i = 1; i < n; i++)
        {
            yx[i + j * d] += alpha * xx[i + j * d];
        }
    }
}

void write_gridfunc2d(pgridfunc2d gf, const char *filename)
{
    FILE *file;
    uint i, j, d;
    double h;
    field *x;

    assert(gf != NULL);

    x = gf->x;
    d = gf->d;
    h = 1.0 / (d - 1.0);

    file = fopen(filename, "w");

    for (i = 0; i < d; i++)
    {
        for (j = 0; j < d; j++)
        {
            fprintf(file, "%.3f\t%.3f\t%+05.3f\n", i * h, j * h, x[i + j * d]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

void init_sine_gridfunc2d(pgridfunc2d gf)
{
    uint i, j, d, n;
    real omega, h, amp;
    field *x;

    assert(gf != NULL);

    x = gf->x;
    d = gf->d;
    n = gf->n;
    h = 1.0 / (n + 1.0);
    omega = 2.0;
    amp = 1.0;

    // Clear left border
    for (i = 0; i < d; i++)
    {
        x[i] = 0.0;
    }

    for (j = 1; j <= n; j++)
    {
        // Clear upper border
        x[j * d] = 0.0;
        for (i = 1; i <= n; i++)
        {
            x[i + j * d] = amp * sin(2.0 * M_PI * omega * (i - 1) * h) *
                           sin(2.0 * M_PI * omega * (j - 1) * h);
        }
        // Clear lower border
        assert(i == n + 1);
        x[i + j * d] = 0.0;
    }

    //Clear right border
    assert(j == n + 1);
    for (i = 0; i < d; i++)
    {
        x[i + j * d] = 0.0;
    }
}

void addeval_5point_stencil_gridfunc2d(pgridfunc2d y, field alpha, pgridfunc2d x)
{
    uint n = y->n;
    uint d = y->d;
    field *xx = x->x;
    field *yy = y->x;

    uint i, j;
    real h;
    field beta;

    assert(x->n == n);
    assert(x->d == d);

    h = 1.0 / (n + 1.0);
    beta = alpha / (h * h);

    for (j = 1; j <= n; j++)
    {
        for (i = 1; i <= n; i++)
        {
            yy[i + j * d] += beta * (4.0 * xx[i + j * d] - xx[(i + 1) + j * d] - xx[(i - 1) + j * d] - xx[i + (j + 1) * d] - xx[i + (j - 1) * d]);
        }
    }
}

#define EPS_SOLVE 1.0e-12
#define MAX_STEPS 1000

void solve_gridfunc2d(pgridfunc2d x)
{
    uint n = x->n;
    uint d = x->d;

    uint i, n2;
    pgridfunc2d r, p, a;
    real norm, error;
    field gamma, lambda, mu;

    n2 = x->d * x->d;
    r = new_gridfunc2d(n);
    p = new_gridfunc2d(n);
    a = new_gridfunc2d(n);

    norm = nrm2(n2, x->x, 1); // Norm of rhs

    // Init residual vector 'r'
    clear_gridfunc2d(r);
    axpy(n2, 1.0, x->x, 1, r->x, 1);               // copy x --> r
    addeval_5point_stencil_gridfunc2d(r, -1.0, x); // r <-- b-Ax

    // Init search direction 'p'
    clear_gridfunc2d(p);
    axpy(n2, 1.0, r->x, 1, p->x, 1); // copy r --> p

    // Compute error
    error = nrm2(n2, r->x, 1);

    i = 0;
    while (error > EPS_SOLVE * norm && i + 1 != MAX_STEPS)
    {
        clear_gridfunc2d(a);
        addeval_5point_stencil_gridfunc2d(a, 1.0, p); // a = A p

        gamma = dot(n2, p->x, 1, a->x, 1); // lambda = <p, r> / <p, a>
        lambda = dot(n2, p->x, 1, r->x, 1) / gamma;

        axpy(n2, lambda, p->x, 1, x->x, 1);  // x = x + lambda p
        axpy(n2, -lambda, a->x, 1, r->x, 1); // r = r - lambda a

        mu = dot(n2, r->x, 1, a->x, 1) / gamma; // p = r - mu p
        scal(n2, -mu, p->x, 1);
        axpy(n2, 1.0, r->x, 1, p->x, 1);

        // Compute error
        error = nrm2(n2, r->x, 1);

        i++;
    }

    // printf("%d steps needed\n\n", i);

    del_gridfunc2d(r);
    del_gridfunc2d(p);
    del_gridfunc2d(a);
}

real diffnorm_gridfunc2d(pgridfunc2d a, pgridfunc2d b)
{
    uint n = a->n;
    uint d = a->d;
    field *ax = a->x;
    field *bx = b->x;

    uint i, j;

    real norm, diff;

    assert(b->n == n);
    assert(b->d == d);

    norm = 0.0;
    for (j = 1; j <= n; j++)
    {
        for (i = 1; i <= n; i++)
        {
            diff = ax[i + j * d] - bx[i + j * d];
            norm += diff * diff;
        }
    }

    return sqrt(norm);
}