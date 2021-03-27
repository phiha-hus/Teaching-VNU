
#include "basic.h"

real c_gamma = 1.0;

struct _gravitation
{
  int planets; // Number of plantes

  real *mass; // array of masses of each planet
  real *y1;   // first component of the position for each planet
  real *y2;   // second component of the position for each planet
};

typedef struct _gravitation gravitation;
typedef gravitation *pgravitation;

pgravitation new_gravitation()
{
  pgravitation grav;

  grav = (pgravitation)allocmem(sizeof(gravitation));

  grav->planets = 2;
  grav->mass = (real *)allocmem(sizeof(real) * grav->planets);
  grav->y1 = (real *)allocmem(sizeof(real) * grav->planets);
  grav->y2 = (real *)allocmem(sizeof(real) * grav->planets);

  // Define masses and position for both planets
  grav->mass[0] = 10.0;
  grav->y1[0] = -0.5;
  grav->y2[0] = -0.2;

  grav->mass[1] = 2.0;
  grav->y1[1] = 0.5;
  grav->y2[1] = 0.5;

  return grav;
}

void del_gravitation(pgravitation grav)
{
  assert(grav != NULL);

  freemem(grav->mass);
  freemem(grav->y1);
  freemem(grav->y2);

  freemem(grav);
}

/**
 * Performs one step of the newton-iteration
 * 
 * 'grav' contains information about the underlying gravitational problem.
 * 'x1', 'x2' denote the current itertion, namely the position of the lagrange point
 * 'f1', 'f2' denote the current force vector at position (x1, x2).
 */
void newton_step(pgravitation grav, real *x1, real *x2,
                 real *f1, real *f2)
{
  real Df[2][2];
  real f[2];
  real diff1, diff2, dist, dist2, dist3, dist5;
  real det;
  real p1, p2;
  uint i;

  Df[0][0] = Df[0][1] = Df[1][0] = Df[1][1] = 0.0;
  f[0] = f[1] = 0;

  for (i = 0; i < grav->planets; i++)
  {
    // compute diff vector for
    diff1 = grav->y1[i] - *x1;
    diff2 = grav->y2[i] - *x2;
    // compute norm squared
    dist2 = diff1 * diff1 + diff2 * diff2;
    dist = sqrt(dist2);
    // compute powers of 3 and 5 of reciprocal norm
    dist3 = (dist * dist2);
    dist5 = 1.0 / (dist3 * dist2);
    dist3 = 1.0 / dist3;

    // evaluate force
    f[0] += grav->mass[i] * diff1 * dist3;
    f[1] += grav->mass[i] * diff2 * dist3;

    // compute jacobian
    Df[0][0] += grav->mass[i] * (3.0 * diff1 * diff1 * dist5 - dist3);
    Df[1][0] += grav->mass[i] * 3.0 * diff2 * diff1 * dist5;
    Df[0][1] += grav->mass[i] * 3.0 * diff1 * diff2 * dist5;
    Df[1][1] += grav->mass[i] * (3.0 * diff2 * diff2 * dist5 - dist3);
  }

  f[0] *= c_gamma;
  f[1] *= c_gamma;
  Df[0][0] *= c_gamma;
  Df[0][1] *= c_gamma;
  Df[1][0] *= c_gamma;
  Df[1][1] *= c_gamma;

  // Write force vector to pointers
  *f1 = f[0];
  *f2 = f[1];

  // compute solution of linear system by Cramer's rule
  det = 1.0 / (Df[0][0] * Df[1][1] - Df[0][1] * Df[1][0]);
  p1 = -(Df[1][1] * f[0] - Df[0][1] * f[1]) * det;
  p2 = -(Df[0][0] * f[1] - Df[1][0] * f[0]) * det;

  //update for the next iteration:
  *x1 += p1;
  *x2 += p2;
}

int main(int argc, char **argv)
{
  pgravitation grav;
  real x1, x2, f1, f2;
  uint i;
  uint steps;

  // number of steps
  steps = 10;

  // init two planetes
  grav = new_gravitation();

  // initial guess of lagrange points
  x1 = 0.0;
  x2 = 0.0;

  for (i = 0; i < steps; i++)
  {
    newton_step(grav, &x1, &x2, &f1, &f2);

    printf("%.15e %.15e     %.3e %.3e\n",
           x1, x2, f1, f2);
  }

  printf("%.15e %.15e\n",
         x1, x2);

  del_gravitation(grav);

  return 0;
}
