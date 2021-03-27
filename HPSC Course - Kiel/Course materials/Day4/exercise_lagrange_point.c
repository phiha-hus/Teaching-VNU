
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
  /**
  * TODO:
  */
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
