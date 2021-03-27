#include "gridfunc1d.h"

pgridfunc1d new_gridfunc1d(uint n)
{
  pgridfunc1d gf;

  gf = (pgridfunc1d)allocmem(sizeof(gridfunc1d));
  gf->n = n;
  gf->d = n + 2;
  gf->x = (field *)allocmem(gf->d * sizeof(field));

  return gf;
}

void del_gridfunc1d(pgridfunc1d gf)
{
  assert(gf != NULL);

  freemem(gf->x);
  freemem(gf);
}

void clear_gridfunc1d(pgridfunc1d gf)
{
  field *x;
  uint d;
  uint i;

  assert(gf != NULL);
  x = gf->x;
  d = gf->d;

  for (i = 0; i < d; i++)
  {
    x[i] = 0.0;
  }
}

void add_gridfunc1d(pgridfunc1d y, field alpha, pgridfunc1d x)
{
  field *xx, *yx;
  uint d, n;
  uint i;

  assert(y != NULL);
  assert(x != NULL);

  xx = x->x;
  yx = y->x;
  n = y->n;
  d = y->d;

  assert(x->d == d);
  assert(x->n == n);

  for (i = 1; i < d; i++)
  {
    yx[i] += alpha * xx[i];
  }
}

void write_gridfunc1d(pgridfunc1d gf, const char *filename)
{
  FILE *file;
  uint i, d;
  double h;
  field *x;

  assert(gf != NULL);

  x = gf->x;
  d = gf->d;
  h = 1.0 / (d - 1.0);

  file = fopen(filename, "w");

  for (i = 0; i < d; i++)
  {
    fprintf(file, "%.3f\t%+05.3f\n", i * h, x[i]);
    fprintf(file, "\n");
  }

  fclose(file);
}

void init_sine_gridfunc1d(pgridfunc1d gf)
{
  uint i, d, n;
  real omega, h, amp;
  field *x;

  assert(gf != NULL);

  x = gf->x;
  d = gf->d;
  n = gf->n;
  h = 1.0 / (n + 1.0);
  omega = 2.0;
  amp = 1.0;

  // Clear left border
  x[0] = 0.0;

  for (i = 1; i < d; i++)
  {
    x[i] = amp * sin(2.0 * M_PI * omega * (i - 1) * h);
  }

  //Clear right border
  x[n + 1] = 0.0;
}

void addeval_3point_stencil_gridfunc1d(pgridfunc1d y, field alpha, pgridfunc1d x)
{
  uint n = y->n;
  uint d = y->d;
  field *xx = x->x;
  field *yy = y->x;

  uint i;
  real h;
  field beta;

  assert(x->n == n);
  assert(x->d == d);

  h = 1.0 / (n + 1.0);
  beta = alpha / (h * h);

  for (i = 1; i <= n; i++)
  {
    yy[i] += beta * (2.0 * xx[i] - xx[i + 1] - xx[i - 1]);
  }
}

#define EPS_SOLVE 1.0e-12
#define MAX_STEPS 1000

void cg_gridfunc1d(pgridfunc1d x)
{
  uint n = x->n;
  uint d = x->d;
  uint i;
  pgridfunc1d r, p, a;
  real norm, error;
  field gamma, lambda, mu;

  r = new_gridfunc1d(n);
  p = new_gridfunc1d(n);
  a = new_gridfunc1d(n);

  norm = nrm2(d, x->x, 1); // Norm of rhs

  // Init residual vector 'r'
  clear_gridfunc1d(r);
  axpy(d, 1.0, x->x, 1, r->x, 1);                // copy x --> r
  addeval_3point_stencil_gridfunc1d(r, -1.0, x); // r <-- b-Ax

  // Init search direction 'p'
  clear_gridfunc1d(p);
  axpy(d, 1.0, r->x, 1, p->x, 1); // copy r --> p

  // Compute error
  error = nrm2(d, r->x, 1);

  i = 0;
  while (error > EPS_SOLVE * norm && i + 1 != MAX_STEPS)
  {
    clear_gridfunc1d(a);
    addeval_3point_stencil_gridfunc1d(a, 1.0, p); // a = A p

    gamma = dot(d, p->x, 1, a->x, 1); // lambda = <p, r> / <p, a>
    lambda = dot(d, p->x, 1, r->x, 1) / gamma;

    axpy(d, lambda, p->x, 1, x->x, 1);  // x = x + lambda p
    axpy(d, -lambda, a->x, 1, r->x, 1); // r = r - lambda a

    mu = dot(d, r->x, 1, a->x, 1) / gamma; // p = r - mu p
    scal(d, -mu, p->x, 1);
    axpy(d, 1.0, r->x, 1, p->x, 1);

    // Compute error
    error = nrm2(d, r->x, 1);

    i++;
  }

  // printf("%d steps needed\n\n", i);

  del_gridfunc1d(r);
  del_gridfunc1d(p);
  del_gridfunc1d(a);
}

void rich_step_1d(real theta, pgridfunc1d b, pgridfunc1d x, pgridfunc1d d)
{
  /**
   *  TODO:
   */
}

void jac_step_1d(real theta, pgridfunc1d b, pgridfunc1d x, pgridfunc1d d)
{
  /**
   *  TODO:
   */
}

void gs_step_1d(pgridfunc1d b, pgridfunc1d x)
{
  /**
   *  TODO:
   */
}

void prolong_gridfunc1d(pgridfunc1d c, pgridfunc1d f)
{
  /**
   *  TODO:
   */
}

void restrict_gridfunc1d(pgridfunc1d f, pgridfunc1d c)
{
  /**
   *  TODO:
   */
}

void mg_gridfunc1d(uint L, pgridfunc1d *x, pgridfunc1d *b, pgridfunc1d *d)
{
  /**
   *  TODO:
   */
}

real diffnorm_gridfunc1d(pgridfunc1d a, pgridfunc1d b)
{
  uint n = a->n;
  uint d = a->d;
  field *ax = a->x;
  field *bx = b->x;

  uint i;

  real norm, diff;

  assert(b->n == n);
  assert(b->d == d);

  norm = 0.0;

  for (i = 1; i < d; i++)
  {
    diff = ax[i] - bx[i];
    norm += diff * diff;
  }

  return sqrt(norm);
}
