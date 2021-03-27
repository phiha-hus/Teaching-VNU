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
    uint d, d2;
    uint i;

    assert(gf != NULL);
    x = gf->x;
    d = gf->d;
    d2 = d * d;

    for (i = 0; i < d2; i++)
    {
        x[i] = 0.0;
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

field dot_gridfunc2d(pgridfunc2d x, pgridfunc2d y)
{
    field *xx, *yx;
    uint d, d2;
    uint i, j;

    field res;

    assert(y != NULL);
    assert(x != NULL);

    xx = x->x;
    yx = y->x;
    d = y->d;
    d2 = d * d;

    assert(x->d == d);

    res = 0.0;

    for (i = 0; i < d2; i++)
    {
        res += xx[i] * yx[i];
    }

    return res;
}

real norm2_gridfunc2d(pgridfunc2d x)
{
    field *xx;
    uint d;

    field res;

    assert(x != NULL);

    xx = x->x;
    d = x->d;

    res = nrm2(d * d, xx, 1);

    return res;
}

void scale_gridfunc2d(pgridfunc2d x, field alpha)
{
    field *xx;
    uint d;

    assert(x != NULL);

    xx = x->x;
    d = x->d;

    scal(d * d, alpha, xx, 1);
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

    omega = 3.0;
    amp = 1.0;

    // Clear left border
    for (i = 0; i < d; i++)
    {
        x[i] = 0.0;
    }

    for (j = 1; j <= n / 4; j++)
    {
        // Clear upper border
        x[j * d] = 0.0;
        for (i = 1; i <= n / 4; i++)
        {
            x[i + j * d] = amp * sin(2.0 * M_PI * omega * (i - 1) * h) *
                           sin(2.0 * M_PI * omega * (j - 1) * h);
        }
        for (; i < d; i++)
        {
            x[i + j * d] = 0.0;
        }
    }

    for (; j < d; j++)
    {
        for (i = 0; i < d; i++)
        {
            x[i + j * d] = 0.0;
        }
    }
}

void addeval_5point_stencil_gridfunc2d(pgridfunc2d y, field alpha, pgridfunc2d x, void *data)
{
    uint n = y->n;
    uint d = y->d;
    field *xx = x->x;
    field *yy = y->x;

    uint i, j;
    real h;
    field beta;

    (void)data;

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

void addeval_cn_wave2d(pgridfunc2d y, field alpha, pgridfunc2d x, void *data)
{
    field c = *((field *)data);

    add_gridfunc2d(y, alpha, x);
    addeval_5point_stencil_gridfunc2d(y, c * alpha, x, NULL);
}

#define EPS_SOLVE 1.0e-12
#define MAX_STEPS 1000

void solve_cg_gridfunc2d(addeval_t A, pgridfunc2d x, void *data)
{
    uint n = x->n;

    uint i;
    pgridfunc2d r, p, a;
    real norm, error;
    field gamma, lambda, mu;

    r = new_gridfunc2d(n);
    p = new_gridfunc2d(n);
    a = new_gridfunc2d(n);

    norm = norm2_gridfunc2d(x); // Norm of rhs

    // Init residual vector 'r'
    clear_gridfunc2d(r);
    add_gridfunc2d(r, 1.0, x); // copy x --> r
    A(r, -1.0, x, data);       // r <-- b-Ax

    // Init search direction 'p'
    clear_gridfunc2d(p);
    add_gridfunc2d(p, 1.0, r); // copy r --> p

    // Compute error
    error = norm2_gridfunc2d(r);

    i = 0;
    while (error > EPS_SOLVE * norm && i + 1 != MAX_STEPS)
    {
        clear_gridfunc2d(a);
        A(a, 1.0, p, data); // a = A p

        gamma = dot_gridfunc2d(p, a); // lambda = <p, r> / <p, a>
        lambda = dot_gridfunc2d(p, r) / gamma;

        add_gridfunc2d(x, lambda, p);  // x = x + lambda p
        add_gridfunc2d(r, -lambda, a); // r = r - lambda a

        mu = dot_gridfunc2d(r, a) / gamma; // p = r - mu p
        scale_gridfunc2d(p, -mu);
        add_gridfunc2d(p, 1.0, r);

        // Compute error
        error = norm2_gridfunc2d(r);

        i++;
    }

    del_gridfunc2d(r);
    del_gridfunc2d(p);
    del_gridfunc2d(a);
}
