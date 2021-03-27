#include "ode.h"
void euler_step(real t, pvector yt, ode_func f, real delta, pvector yt1, void *data)
{
    /*
     * TODO
     */

    uint dim = yt->dim;

    uint i;

    assert(dim == yt1->dim);

    // evaluate function at point 't', current vector 'yt' into 'yt1'
    f(t, yt, yt1, data);

    // update new vector 'yt1'
    for (i = 0; i < dim; i++)
    {
        yt->x[i] += delta * yt1->x[i];
    }    
}

void runge_step(real t, pvector yt, ode_func f, real delta, pvector yt1, void *data)
{
    /*
     * TODO
     */
}

void leapfrog_step(real t, pvector yt, ode_func f, real delta, pvector yt1, void *data)
{
    /*
     * TODO
     */
}

void crank_nicolson_step(real t, pvector yt, ode_func f, real delta, pvector yt1, void *data)
{
    /*
     * TODO
     */
}
