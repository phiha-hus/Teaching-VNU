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
    real gamma;
    pvector yt, yt1;
    FILE *fp;
    char filename[100];

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

    delta = 1.0e-5;
    printf("Computing reference value in endpoint with rk-classic, delta = %.2e\n", delta);

    /**
     *  TODO: 
     */

    printf(
        "Simulating spring-mass system via explicit Euler method:\n");

    delta = 5.0e-1;
    while (delta >= 1.0e-4)
    {

        init_mass_spring(yt);

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
        /**
        * TODO:
        */

        delta *= 0.5;
    }

    printf(
        "Simulating spring-mass system via classical Runge-Kutta method:\n");

    /**
     * TODO:
     */

    printf(
        "Simulating spring-mass system via Adams-Bashforth method of order 3:\n");

    /**
     * TODO:
     */

    del_vector(yt);
    del_vector(yt1);

    return 0;
}
