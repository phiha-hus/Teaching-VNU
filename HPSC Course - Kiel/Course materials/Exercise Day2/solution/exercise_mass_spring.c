#include "basic.h"
#include "linalg.h"
#include "ode.h"

/**
 * Define things for mass-spring system
 */

real mass_spring_c = 1.0;
real mass_spring_m = 0.5;

/**
 * Setup some initial values for the mass-spring system
 */
void init_mass_spring(pvector y0)
{
    assert(y0->dim == 2);

    y0->x[0] = 1.0; // initial positon at right border
    y0->x[1] = 0.0; // initial velocity pointing to the left
}

/**
 * right-hand-side for the mass-spring system
 */
void mass_spring_func(real t, pvector yt, pvector yt1, void *data)
{
    real c_m = *((real *)data); // the constant c/m appearing in the equations

    yt1->x[0] = yt->x[1];
    yt1->x[1] = -c_m * yt->x[0];
}

int main(int argc, char const *argv[])
{

    uint dim;
    real a, b;
    real delta;
    real t;
    uint steps;
    real mass_spring_c_m;
    pvector yt, yt1, yt_ref;
    pvector *dy;
    uint m, i;
    real *w;
    FILE *fp;
    char filename[100];
    real error, last_error;

    printf(
        "\n"
        "##############################################################################\n"
        "#                      Exercise: mass-spring system                          #\n"
        "##############################################################################\n"
        "\n");

    /****************************************************************************
     *  Mass-spring system
     *****************************************************************************/

    dim = 2;
    mass_spring_c_m = mass_spring_c / mass_spring_m;
    a = 0.0;
    b = 16.0;

    // Initialize problem
    yt = new_vector(dim);
    yt1 = new_vector(dim);
    yt_ref = new_vector(dim);

    printf("Computing reference value in endpoint with rk-classic, delta = %.2e\n", delta);

    delta = 5.0e-5;
    init_mass_spring(yt_ref);

    // start simulation from 't = a' to 't = b' with stepwidth 'delta'
    t = a;
    while (t + delta < b)
    {
        rk_classic_step(t, yt_ref, (ode_func)mass_spring_func, delta, yt1, &mass_spring_c_m);
        t += delta;
    }
    rk_classic_step(t, yt_ref, (ode_func)mass_spring_func, b - t, yt1, &mass_spring_c_m);

    printf(
        "Simulating spring-mass system via explicit Euler method:\n");

    delta = 5.0e-1;
    last_error = 0.0;
    while (delta >= 1.0e-4)
    {

        init_mass_spring(yt);
        error = 0.0;

        // start simulation from 't = a' to 't = b' with stepwidth 'delta'
        t = a;
        sprintf(filename, "data/ode_euler_delta%.4f.dat", delta);
        fp = fopen(filename, "w");
        assert(fp != NULL);
        fprintf(fp, "%.4f\t%.5e\t%.5e\n", t, yt->x[0], yt->x[1]);
        printf("  Starting simulation with delta = %.2e\n", delta);
        while (t + delta < b)
        {
            euler_step(t, yt, (ode_func)mass_spring_func, delta, yt1, &mass_spring_c_m);
            t += delta;
            fprintf(fp, "%.4f\t%.5e\t%.5e\n", t, yt->x[0], yt->x[1]);
        }
        euler_step(t, yt, (ode_func)mass_spring_func, b - t, yt1, &mass_spring_c_m);
        t += b - t;
        fprintf(fp, "%.4f\t%.5e\t%.5e\n", t, yt->x[0], yt->x[1]);
        fclose(fp);
        printf("  final position: %.3f, final velocity: %.3f\n", yt->x[0], yt->x[1]);

        // Compute error
        for (i = 0; i < dim; i++)
        {
            error = fabs(yt_ref->x[i] - yt->x[i]) > error ? fabs(yt_ref->x[i] - yt->x[i]) : error;
        }
        printf("  max error: %.5e (x %.2f)\n", error, last_error / error);
        printf("\n");
        last_error = error;

        delta *= 0.5;
    }

    printf(
        "Simulating spring-mass system via classical Runge-Kutta method:\n");

    delta = 5.0e-1;
    last_error = 0.0;
    while (delta >= 1.0e-4)
    {

        init_mass_spring(yt);
        error = 0.0;

        // start simulation from 't = a' to 't = b' with stepwidth 'delta'
        t = a;
        sprintf(filename, "data/ode_rk_classic_delta%.4f.dat", delta);
        fp = fopen(filename, "w");
        assert(fp != NULL);
        fprintf(fp, "%.4f\t%.5e\t%.5e\n", t, yt->x[0], yt->x[1]);
        printf("  Starting simulation with delta = %.2e\n", delta);
        while (t + delta < b)
        {
            rk_classic_step(t, yt, (ode_func)mass_spring_func, delta, yt1, &mass_spring_c_m);
            t += delta;
            fprintf(fp, "%.4f\t%.5e\t%.5e\n", t, yt->x[0], yt->x[1]);
        }
        rk_classic_step(t, yt, (ode_func)mass_spring_func, b - t, yt1, &mass_spring_c_m);
        t += b - t;
        fprintf(fp, "%.4f\t%.5e\t%.5e\n", t, yt->x[0], yt->x[1]);
        fclose(fp);
        printf("  final position: %.3f, final velocity: %.3f\n", yt->x[0], yt->x[1]);

        // Compute error
        for (i = 0; i < dim; i++)
        {
            error = fabs(yt_ref->x[i] - yt->x[i]) > error ? fabs(yt_ref->x[i] - yt->x[i]) : error;
        }
        printf("  max error: %.5e (x %.2f)\n", error, last_error / error);
        printf("\n");
        last_error = error;

        delta *= 0.5;
    }

    printf(
        "Simulating spring-mass system via Adams-Bashforth method of order 3:\n");

    m = 3;
    delta = 5.0e-1;
    last_error = 0.0;
    while (delta >= 1.0e-4)
    {

        // start simulation from 't = a' to 't = b' with stepwidth 'delta'
        t = a;

        init_mass_spring(yt);
        error = 0.0;
        adams_bashforth_init(t, yt, m, &dy, &w, (ode_func)mass_spring_func, delta, &mass_spring_c_m);
        t += (m + 1) * delta;
        i = 0;

        sprintf(filename, "data/ode_adams_bashforth_3_delta%.4f.dat", delta);
        fp = fopen(filename, "w");
        assert(fp != NULL);
        fprintf(fp, "%.4f\t%.5e\t%.5e\n", t, yt->x[0], yt->x[1]);
        printf("  Starting simulation with delta = %.2e\n", delta);
        while (t < b)
        {
            adams_bashforth_step(t, yt, m, dy, w, i, (ode_func)mass_spring_func, delta, &mass_spring_c_m);
            i = (i + 1) % (m + 1);
            t += delta;
            fprintf(fp, "%.4f\t%.5e\t%.5e\n", t, yt->x[0], yt->x[1]);
        }

        fclose(fp);
        printf("  final position: %.3f, final velocity: %.3f\n", yt->x[0], yt->x[1]);

        // Compute error
        for (i = 0; i < dim; i++)
        {
            error = fabs(yt_ref->x[i] - yt->x[i]) > error ? fabs(yt_ref->x[i] - yt->x[i]) : error;
        }
        printf("  max error: %.5e (x %.2f)\n", error, last_error / error);
        printf("\n");
        last_error = error;

        delta *= 0.5;
    }

    del_vector(yt);
    del_vector(yt1);

    return 0;
}
