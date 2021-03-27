
#include "basic.h"
#include "gridfunc2d.h"

void clear_omp_gridfunc2d(pgridfunc2d gf)
{
    field *x;
    uint d, d2;
    uint i;

    assert(gf != NULL);
    x = gf->x;
    d = gf->d;
    d2 = d * d;

    // TODO: parallelize this code using OpenMP for-directive

    for (i = 0; i < d2; i++)
    {
        x[i] = 0.0;
    }
}

void clear_omp_api_gridfunc2d(pgridfunc2d gf)
{
    field *x;
    uint d, d2;
    uint i;

    assert(gf != NULL);
    x = gf->x;
    d = gf->d;
    d2 = d * d;

    // TODO: parallelize this code using OpenMP API-functions for worksharing

    for (i = 0; i < d2; i++)
    {
        x[i] = 0.0;
    }
}

field dot_omp_gridfunc2d(pgridfunc2d x, pgridfunc2d y)
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

    // TODO: parallelize this code using OpenMP reduction clause

    res = 0.0;
    for (i = 0; i < d2; i++)
    {
        res += xx[i] * yx[i];
    }

    return res;
}

field dot_omp_manual_gridfunc2d(pgridfunc2d x, pgridfunc2d y)
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

    // TODO: parallelize this code using local variables for reduction and sychronization

    res = 0.0;
    for (i = 0; i < d2; i++)
    {
        res += xx[i] * yx[i];
    }

    return res;
}

void addeval_omp_5point_stencil_gridfunc2d(pgridfunc2d y, field alpha, pgridfunc2d x, void *data)
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

    // TODO: parallelize this code using OpenMP for-directive

    for (j = 1; j <= n; j++)
    {
        for (i = 1; i <= n; i++)
        {
            yy[i + j * d] += beta * (4.0 * xx[i + j * d] - xx[(i + 1) + j * d] - xx[(i - 1) + j * d] - xx[i + (j + 1) * d] - xx[i + (j - 1) * d]);
        }
    }
}

int main(int argc, char const *argv[])
{

    pgridfunc2d x, y;
    pstopwatch sw;
    uint n;
    uint iter;
    uint i;
    real t, t2;
    field val;

    n = 10000;
    iter = 20;

    omp_set_nested(1);

    sw = new_stopwatch();

    x = new_gridfunc2d(n);
    y = new_gridfunc2d(n);

    /*************************
     * Test dot_gridfunc2d
     *************************/

    printf("Measure %d evaluation of standard clear_gridfunc2d:\n", iter);
    start_stopwatch(sw);

    for (i = 0; i < iter; i++)
    {
        clear_gridfunc2d(y);
    }

    t2 = stop_stopwatch(sw);
    printf("  %.3f s\n", t2);

    printf("Measure %d evaluation of parallel (API) clear_gridfunc2d:\n", iter);
    start_stopwatch(sw);

    for (i = 0; i < iter; i++)
    {
        clear_omp_api_gridfunc2d(y);
    }

    t = stop_stopwatch(sw);
    printf("  %.3f s\n", t);
    printf("Speedup: %.2f\n", t2 / t);

    printf("Measure %d evaluation of parallel (FOR) clear_gridfunc2d:\n", iter);
    start_stopwatch(sw);

    for (i = 0; i < iter; i++)
    {
        clear_omp_gridfunc2d(y);
    }

    t = stop_stopwatch(sw);
    printf("  %.3f s\n", t);
    printf("Speedup: %.2f\n", t2 / t);

    printf("\n##############################################\n\n");

    /*************************
     * Test dot_gridfunc2d
     *************************/

    init_sine_gridfunc2d(x);
    init_sine_gridfunc2d(y);

    printf("Measure %d evaluation of standard dot_gridfunc2d:\n", iter);
    start_stopwatch(sw);

    for (i = 0; i < iter; i++)
    {
        val = dot_gridfunc2d(x, y);
    }

    t2 = stop_stopwatch(sw);
    printf("  %.3f s\n", t2);
    printf("val = %.5e\n", val);

    printf("Measure %d evaluation of parallel (reduce) dot_gridfunc2d:\n", iter);
    start_stopwatch(sw);

    for (i = 0; i < iter; i++)
    {
        val = dot_omp_gridfunc2d(x, y);
    }

    t = stop_stopwatch(sw);
    printf("  %.3f s\n", t);
    printf("val = %.5e\n", val);
    printf("Speedup: %.2f\n", t2 / t);

    printf("Measure %d evaluation of parallel (manual) dot_gridfunc2d:\n", iter);
    start_stopwatch(sw);

    for (i = 0; i < iter; i++)
    {
        val = dot_omp_manual_gridfunc2d(x, y);
    }

    t = stop_stopwatch(sw);
    printf("  %.3f s\n", t);
    printf("val = %.5e\n", val);
    printf("Speedup: %.2f\n", t2 / t);

    printf("\n##############################################\n\n");

    /*************************
     * Test addeval_5point_stencil_gridfunc2d
     *************************/

    init_sine_gridfunc2d(x);
    clear_gridfunc2d(y);

    printf("Measure %d evaluation of standard addeval_5point_stencil:\n", iter);
    start_stopwatch(sw);

    for (i = 0; i < iter; i++)
    {
        addeval_5point_stencil_gridfunc2d(y, 1.0, x, NULL);
    }

    t2 = stop_stopwatch(sw);
    printf("  %.3f s\n", t2);

    printf("Measure %d evaluation of parallel addeval_5point_stencil:\n", iter);
    start_stopwatch(sw);

    for (i = 0; i < iter; i++)
    {
        addeval_omp_5point_stencil_gridfunc2d(y, -1.0, x, NULL);
    }

    t = stop_stopwatch(sw);
    printf("  %.3f s\n", t);
    printf("Speedup: %.2f\n", t2 / t);

    printf("norm of y is %.5e\n", norm2_gridfunc2d(y));

    del_stopwatch(sw);
    del_gridfunc2d(x);
    del_gridfunc2d(y);

    return 0;
}
