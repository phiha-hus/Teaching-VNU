#include "basic.h"
#include "linalg.h"
#include "ode.h"

/**
 * Define things for gravitational N-body problem
 */

/**
 * Setup some initial values for the gravitation problem
 */
void init_gravitation(pvector y0)
{
    uint dim = y0->dim;
    // Vector has length 7*n, 3*n for positions, 3*n for velocities, n for masses
    uint n = dim / 7;
    uint i;
    real scale = 1.0 / RAND_MAX;

    srand(42);

    for (i = 0; i < n; i++)
    {
        // set some random positions
        y0->x[i + 0 * n] = (real)rand() * scale;
        y0->x[i + 1 * n] = (real)rand() * scale;
        y0->x[i + 2 * n] = (real)rand() * scale;
        // set velocities to zero
        y0->x[i + 3 * n] = 0.0;
        y0->x[i + 4 * n] = 0.0;
        y0->x[i + 5 * n] = 0.0;
        // set some random masses
        y0->x[i + 6 * n] = 10.0 + 10.0 * ((real)rand() * scale);
    }
}

/**
 * right-hand-side for the gravitation problemm
 */
void gravitation_func(real t, pvector yt, pvector yt1, void *data)
{
    uint dim = yt->dim;
    // Vector has length 7*n, 3*n for positions, 3*n for velocities, n for masses
    uint n = dim / 7;
    real gamma = *((real *)data);
    uint i, j;
    real sum[3];
    real dist[3];
    real norm;
    real m;

    assert(yt1->dim == dim);

    for (i = 0; i < n; i++)
    {
        // x'(t) = v(t)
        yt1->x[i + 0 * n] = yt->x[i + 3 * n];
        yt1->x[i + 1 * n] = yt->x[i + 4 * n];
        yt1->x[i + 2 * n] = yt->x[i + 5 * n];

        // compute forces
        sum[0] = 0.0;
        sum[1] = 0.0;
        sum[2] = 0.0;

        for (j = 0; j < n; j++)
        {
            if (i != j)
            {
                // Get the mass m_j
                m = yt->x[j + 6 * n];

                // Distance vector
                dist[0] = yt->x[j + 0 * n] - yt->x[i + 0 * n];
                dist[1] = yt->x[j + 1 * n] - yt->x[i + 1 * n];
                dist[2] = yt->x[j + 2 * n] - yt->x[i + 2 * n];

                norm = dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];
                norm = 1.0 / sqrt(norm);
                norm = m * norm * norm * norm;

                sum[0] += dist[0] * norm;
                sum[1] += dist[1] * norm;
                sum[2] += dist[2] * norm;
            }
        }

        // set velocities
        yt1->x[i + 3 * n] = gamma * sum[0];
        yt1->x[i + 4 * n] = gamma * sum[1];
        yt1->x[i + 5 * n] = gamma * sum[2];

        // set "new mass" to zero.
        yt1->x[i + 6 * n] = 0.0;
    }
}

void print_state_gravitation(pvector y)
{
    uint dim = y->dim;
    // Vector has length 7*n, 3*n for positions, 3*n for velocities, n for masses
    uint n = dim / 7;
    uint i;

    printf("%d bodies\n", n);
    printf("                  mass      |                    position                   |                    velocity\n");
    for (i = 0; i < n && i < 10; i++)
    {
        printf("mass %2d:\t%.3e\t(%+.5e, %+.5e, %+.5e)\t(%+.5e, %+.5e, %+.5e)\n", i, y->x[i + 6 * n],
               y->x[i + 0 * n], y->x[i + 1 * n], y->x[i + 2 * n],
               y->x[i + 3 * n], y->x[i + 4 * n], y->x[i + 5 * n]);
    }
    printf("\n");
}

int main(int argc, char const *argv[])
{

    uint dim;
    real a, b;
    real delta;
    real t;
    uint steps;
    real gamma;
    pvector yt, yt1;
    FILE *fp;
    char filename[100];

    printf(
        "\n"
        "##############################################################################\n"
        "#                               Exercise 02                                  #\n"
        "##############################################################################\n"
        "\n");

    /****************************************************************************
     *  Many body problem
     *****************************************************************************/

    printf(
        "Simulating gravitational N-Body Problem via explicit Euler method:\n");

    dim = 20;
    gamma = 1.0e-5;
    a = 0.0;
    b = 2.0;

    // Initialize problem
    yt = new_vector(dim * 7);
    yt1 = new_vector(dim * 7);

    delta = 1.0e-1;
    while (delta >= 1.0e-5)
    {

        init_gravitation(yt);
        printf("Initial configuration (delta = %.2e):\n", delta);
        print_state_gravitation(yt);

        // start simulation from 't = a' to 't = b' with stepwidth 'delta'
        t = a;
        while (t < b)
        {
            euler_step(t, yt, (ode_func)gravitation_func, delta, yt1, &gamma);
            t += delta;
        }
        printf("Final configuration (delta = %.2e):\n", delta);
        print_state_gravitation(yt);
        printf("\n\n");

        delta *= 0.1;
    }

    del_vector(yt);
    del_vector(yt1);

    return 0;
}
