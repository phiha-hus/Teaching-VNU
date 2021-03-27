
#include "basic.h"
#include "gridfunc2d.h"

real wave_c = 0.5;

void explicit_euler(pgridfunc2d x, pgridfunc2d v, real delta, real t_start, real t_end)
{
    char filename[255];
    real t;
    uint i;

    printf("Explicit Euler:\n");

    t = t_start;

    i = 0;
    while (t < t_end)
    {
        // print out current time 't'
        printf("\rt = %.3f s", t);

        // TODO: implement the Euler step for the wave equation.

        if (!(i % 100))
        {
            sprintf(filename, "data/wave2d_euler_%05d.dat", i);
            write_gridfunc2d(x, filename);
        };

        t += delta;
        i++;
    }
    printf("\n");
}

void leapfrog(pgridfunc2d x, pgridfunc2d v, real delta, real t_start, real t_end)
{
    char filename[255];
    real t;
    uint i;

    printf("Leapfrog:\n");

    t = t_start;

    // Initialization phase
    i = 0;
    // print out current time 't'
    printf("\rt = %.3f s", t);

    // TODO: Initialization of leapfrog scheme

    sprintf(filename, "data/wave2d_leapfrog_%05d.dat", i);
    write_gridfunc2d(x, filename);
    t += delta;
    i++;

    while (t < t_end)
    {
        // print out current time 't'
        printf("\rt = %.3f s", t);

        // TODO: implement the Leapfrog step for the wave equation.

        if (!(i % 100))
        {
            sprintf(filename, "data/wave2d_leapfrog_%05d.dat", i);
            write_gridfunc2d(x, filename);
        };

        t += 2.0 * delta;
        i++;
    }
    printf("\n");
}

int main(int argc, char const *argv[])
{
    pgridfunc2d x, v;
    uint n;
    real delta;
    real t_start;
    real t_end;

    // number of Grid points per dimension
    n = 50;
    // timestep size
    delta = 1.0 / (n * n);
    // time interval
    t_start = 0.0;
    t_end = 4.0;

    x = new_gridfunc2d(n);
    v = new_gridfunc2d(n);

    // Intial values for x at t=t_0, v at t=t_0 is zero.
    init_sine_gridfunc2d(x);
    clear_gridfunc2d(v);

    // perform computations with explicit euler
    explicit_euler(x, v, delta, t_start, t_end);

    // Intial values for x at t=t_0, v at t=t_0 is zero.
    init_sine_gridfunc2d(x);
    clear_gridfunc2d(v);

    // perform computations with leapfrog
    // leapfrog needs half the timestep size, because
    // it makes the 'leaps'.
    leapfrog(x, v, 0.5 * delta, t_start, t_end);

    // cleanup
    del_gridfunc2d(x);
    del_gridfunc2d(v);

    return 0;
}
