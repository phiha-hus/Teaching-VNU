#include "ode.h"
void euler_step(real t, pvector yt, ode_func f, real delta, pvector yt1, void *data)
{
    uint dim = yt->dim;

    uint i;

    assert(dim == yt1->dim);

    // evaluate function at point 't', current vector 'yt' into 'yt1'
    f(t, yt, yt1, data);

    // update new vector 'yt1'
    for (i = 0; i < dim; i++)
    {
        yt->x[i] = yt->x[i] + delta * yt1->x[i];
    }
}

void runge_step(real t, pvector yt, ode_func f, real delta, pvector yt1, void *data)
{
    uint dim = yt->dim;

    pvector ym;
    uint i;

    assert(dim == yt1->dim);

    ym = new_vector(dim);

    // evaluate function at point 't', current vector 'yt' into 'yt1'
    f(t, yt, yt1, data);

    // update new vector 'ym'
    for (i = 0; i < dim; i++)
    {
        ym->x[i] = yt->x[i] + 0.5 * delta * yt1->x[i];
    }

    // evaluate function at point 't+delta/2', current vector 'ym' into 'yt1'
    f(t, ym, yt1, data);

    // update new vector 'yt'
    for (i = 0; i < dim; i++)
    {
        yt->x[i] = yt->x[i] + delta * yt1->x[i];
    }

    del_vector(ym);
}

void leapfrog_step(real t, pvector yt, ode_func f, real delta, pvector yt1, void *data)
{
    uint dim = yt->dim;
    real c_mass_spring = *((real *)data);

    (void)f;
    (void)yt1;

    assert(dim == yt1->dim);
    assert(dim == 2);

    yt->x[0] = yt->x[0] + 2.0 * delta * yt->x[1];
    yt->x[1] = yt->x[1] - 2.0 * delta * c_mass_spring * yt->x[0];
}

void crank_nicolson_step(real t, pvector yt, ode_func f, real delta, pvector yt1, void *data)
{
    uint dim = yt->dim;
    real c_mass_spring = *((real *)data);

    uint i;
    real onepdelta, onemdelta;

    (void)f;

    assert(dim == yt1->dim);
    assert(dim == 2);

    onepdelta = 1.0 + 0.25 * delta * delta * c_mass_spring;
    onemdelta = 1.0 - 0.25 * delta * delta * c_mass_spring;

    yt1->x[0] = onemdelta * yt->x[0] + delta * yt->x[1];
    yt1->x[1] = onemdelta * yt->x[1] - delta * c_mass_spring * yt->x[0];

    yt->x[0] = yt1->x[0] / onepdelta;
    yt->x[1] = yt1->x[1] / onepdelta;
}

void rk_classic_step(real t, pvector yt, ode_func f, real delta, pvector yt1, void *data)
{
    uint dim = yt->dim;

    pvector tmp, k1, k2, k3, k4;

    tmp = new_vector(dim);
    k1 = new_vector(dim);
    k2 = new_vector(dim);
    k3 = new_vector(dim);
    k4 = new_vector(dim);

    // Evaluate intermediate results k1 ... k4

    //k1
    f(t, yt, k1, data);

    //k2
    clear_vector(tmp);
    add_vector(tmp, 1.0, yt);
    add_vector(tmp, 0.5 * delta, k1);
    f(t + 0.5 * delta, tmp, k2, data);

    //k3
    clear_vector(tmp);
    add_vector(tmp, 1.0, yt);
    add_vector(tmp, 0.5 * delta, k2);
    f(t + 0.5 * delta, tmp, k3, data);

    //k4
    clear_vector(tmp);
    add_vector(tmp, 1.0, yt);
    add_vector(tmp, delta, k3);
    f(t + delta, tmp, k4, data);

    //combine final result
    add_vector(yt, delta / 6.0, k1);
    add_vector(yt, delta / 3.0, k2);
    add_vector(yt, delta / 3.0, k3);
    add_vector(yt, delta / 6.0, k4);

    del_vector(tmp);
    del_vector(k1);
    del_vector(k2);
    del_vector(k3);
    del_vector(k4);
}

/**
 * ATTENTION: SOME ERROR IN THIS METHOD, NOT THE CORRECT RATE OF CONVERGENCE
 */

void adams_bashforth_init(real t, pvector yt, uint m, pvector **dy, real **w, ode_func f, real delta, void *data)
{
    uint dim = yt->dim;

    pvector yt1;
    uint i;

    assert(m == 3);

    yt1 = new_vector(dim);
    *dy = (pvector *)allocmem((m + 1) * sizeof(pvector));
    *w = (real *)allocmem((m + 1) * sizeof(real));

    for (i = 0; i <= m; i++)
    {
        (*dy)[i] = new_vector(dim);
        f(t, yt, (*dy)[i], data);
        rk_classic_step(t, yt, f, delta, yt1, data);
        t += delta;
    }

    (*w)[0] = 55.0 / 24.0;
    (*w)[1] = -59.0 / 24.0;
    (*w)[2] = 37.0 / 24.0;
    (*w)[3] = -9.0 / 24.0;

    del_vector(yt1);
}

void adams_bashforth_step(real t, pvector yt, uint m, pvector *dy, real *w, uint i, ode_func f, real delta, void *data)
{
    assert(m > 0);
    assert(i <= m);

    uint dim = yt->dim;

    uint j, k;

    for (j = 0; j <= m; j++)
    {
        // get index for old values of derivatives
        k = (i + (m + 1) - j) % (m + 1);
        add_vector(yt, delta * w[j], dy[k]);
    }

    k = (i + 1) % (m + 1);
    f(t + delta, yt, dy[k], data);
}