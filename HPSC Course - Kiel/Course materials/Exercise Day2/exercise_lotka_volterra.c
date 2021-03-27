#include "basic.h"
#include "linalg.h"
#include "ode.h"

void init_lotka_volterra(pvector y)
{
    uint dim = y->dim;

    assert(dim == 2);

    y->x[0] = 0.5;
    y->x[1] = 0.5;
}

void lotka_volterra_func(real t, pvector yt, pvector yt1, void *data)
{
    real *params = (real *)data;

    /**
     * TODO:
     */
}

int main(int argc, char const *argv[])
{

    uint dim;
    real a, b;
    real delta;
    real params_lotka_volterra[4];
    real t;
    pvector yt, yt1;
    FILE *fp;
    char filename[100];

    printf(
        "\n"
        "##############################################################################\n"
        "#                        Exercise Lotka-Volterra                             #\n"
        "##############################################################################\n"
        "\n");

    dim = 2;
    params_lotka_volterra[0] = 0.2;  // alpha;
    params_lotka_volterra[1] = 0.2;  // beta;
    params_lotka_volterra[2] = 0.25; // gamma;
    params_lotka_volterra[3] = 0.15; // delta;
    a = 0.0;
    b = 250.0;

    // Initialize problem
    yt = new_vector(dim);
    yt1 = new_vector(dim);

    printf(
        "Simulating Lotka-Volterra via explicit Euler method:\n");

    delta = 4.0e-1;
    while (delta >= 1.0e-4)
    {

        init_lotka_volterra(yt);

        // start simulation from 't = a' to 't = b' with stepwidth 'delta'
        t = a;
        sprintf(filename, "data/LV_euler_delta%.4f.dat", delta);
        fp = fopen(filename, "w");
        assert(fp != NULL);
        fprintf(fp, "%.4f\t%.5e\t%.5e\n", t, yt->x[0], yt->x[1]);
        while (t < b)
        {
            euler_step(t, yt, (ode_func)lotka_volterra_func, delta, yt1, (void *)params_lotka_volterra);
            t += delta;
            fprintf(fp, "%.4f\t%.5e\t%.5e\n", t, yt->x[0], yt->x[1]);
        }
        fclose(fp);

        delta *= 0.25;
    }

    printf(
        "Simulating Lotka-Volterra via Runge's method:\n");

    delta = 4.0e-1;
    while (delta >= 1.0e-4)
    {

        init_lotka_volterra(yt);

        // start simulation from 't = a' to 't = b' with stepwidth 'delta'
        t = a;

        sprintf(filename, "data/LV_runge_delta%.4f.dat", delta);
        fp = fopen(filename, "w");
        assert(fp != NULL);
        fprintf(fp, "%.4f\t%.5e\t%.5e\n", t, yt->x[0], yt->x[1]);
        while (t < b)
        {
            runge_step(t, yt, (ode_func)lotka_volterra_func, delta, yt1, (void *)params_lotka_volterra);
            t += delta;
            fprintf(fp, "%.4f\t%.5e\t%.5e\n", t, yt->x[0], yt->x[1]);
        }
        fclose(fp);

        delta *= 0.25;
    }

    printf(
        "Simulating Lotka-Volterra via classical RK4 method:\n");

    delta = 4.0e-1;
    while (delta >= 1.0e-4)
    {

        init_lotka_volterra(yt);

        // start simulation from 't = a' to 't = b' with stepwidth 'delta'
        t = a;

        sprintf(filename, "data/LV_rk4_delta%.4f.dat", delta);
        fp = fopen(filename, "w");
        assert(fp != NULL);
        fprintf(fp, "%.4f\t%.5e\t%.5e\n", t, yt->x[0], yt->x[1]);
        while (t < b)
        {
            rk_classic_step(t, yt, (ode_func)lotka_volterra_func, delta, yt1, (void *)params_lotka_volterra);
            t += delta;
            fprintf(fp, "%.4f\t%.5e\t%.5e\n", t, yt->x[0], yt->x[1]);
        }
        fclose(fp);

        delta *= 0.25;
    }

    del_vector(yt);
    del_vector(yt1);

    return 0;
}
