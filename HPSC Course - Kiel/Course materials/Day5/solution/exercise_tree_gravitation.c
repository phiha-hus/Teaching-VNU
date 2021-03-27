

#include "basic.h"
#include "cluster.h"

/* ------------------------------------------------------------
 * Evaluate gravitation potential
 * ------------------------------------------------------------ */

/**
 * Evaluate the potential between points (x1, x2, x3) and (y1, y2, y3)
 */
static real
potential(real x1, real x2, real x3,
          real y1, real y2, real y3)
{
  real norm2, norm;

  norm2 = sqr(y1 - x1) + sqr(y2 - x2) + sqr(y3 - x3);
  norm = sqrt(norm2);

  return 1.0 / norm;
}

/**
 * Eval the gravitational potential at point (x1,x2,x3) for
 * one cluster 's' using the interpolation approximation.
 */
real eval_approximation(real x1, real x2, real x3, pcluster s)
{
  int m = s->m;               // Order of interpolation
  const real *xi1 = s->xi[0]; // Interpolation points in x-direction
  const real *xi2 = s->xi[1]; // Interpolation points in y-direction
  const real *xi3 = s->xi[2]; // Interpolation points in z-direction
  uint nu1, nu2, nu3, i1, i2, i3;

  real result;

  result = 0.0;
  for (nu1 = 0; nu1 <= m; nu1++)
  {
    i1 = nu1;

    for (nu2 = 0; nu2 <= m; nu2++)
    {
      i2 = i1 * (m + 1) + nu2;

      for (nu3 = 0; nu3 <= m; nu3++)
      {
        i3 = i2 * (m + 1) + nu3;

        result += potential(x1, x2, x3, xi1[nu1], xi2[nu2], xi3[nu3]) * s->z[i3];
      }
    }
  }

  return result;
}

/**
 * Eval the gravitational potential at point (x1,x2,x3) for
 * a single cluster 's' with the direct evaluation.
 * The planets are permuted by 's->idx'
 */
real eval_full(real x1, real x2, real x3, pcluster s, const state *st)
{
  const uint size = s->size;
  const int *idx = s->idx;
  const real *y1 = st->y[0];
  const real *y2 = st->y[1];
  const real *y3 = st->y[2];
  uint i;

  real result;

  result = 0.0;
  for (i = 0; i < size; i++)
  {
    result += potential(x1, x2, x3, y1[idx[i]], y2[idx[i]], y3[idx[i]]);
  }

  return result;
}

/**
 * Eval the gravitational potential at point (x1,x2,x3) for
 * all planets 1, ..., st->n using the approxiations.
 */
real eval_cluster(real x1, real x2, real x3,
                  cluster *s, const state *st, real eta)
{
  real dist2, diam2;
  real result;

  uint i;

  // compute the squared distant between x and s.
  dist2 = dist2_cluster(x1, x2, x3, s);
  // compute the squared diameter of 's'
  diam2 = diam2_cluster(s);

  result = 0.0;

  // Admissible case --> evaluate approximation
  if (diam2 <= eta * eta * dist2)
  {
    result += eval_approximation(x1, x2, x3, s);
  }
  // Cluster is inadmissible, but has sons --> recursion
  else if (s->sons > 0)
  {
    for (i = 0; i < s->sons; i++)
      result += eval_cluster(x1, x2, x3, s->son[i], st, eta);
  }
  // cluster is inadmissible and has _NO_ sons --> direct evaluation
  else
  {
    result += eval_full(x1, x2, x3, s, st);
  }

  return result;
}

/**
 * Eval the gravitational potential at point (x1,x2,x3) for
 * all planets 1, ... st->n with the direct evaluation.
 */
real eval_potential(real x1, real x2, real x3,
                    const state *st)
{
  real result;
  int n = st->n;
  const real *y1 = st->y[0];
  const real *y2 = st->y[1];
  const real *y3 = st->y[2];
  int i;

  result = 0.0;

  for (i = 0; i < n; i++)
    result += potential(x1, x2, x3, y1[i], y2[i], y3[i]);

  return result;
}

int main()
{
  state *st;
  interpolation *in;
  cluster *s;
  real bmin[3], bmax[3];
  real *x1, *x2, *x3, *phi1, *phi2;
  real error, maxerror;
  real eta;
  int n, m, n_test;
  int i;

  /* Initialize random number generator */
  srand(42);

  /* Number of planets */
  n = 262144;

  /* Number of test positions */
  n_test = 1000;

  /* Degree of interpolation */
  m = 3;

  /* Admissibility parameter */
  eta = 2.0;

  printf("Evaluation %d times the potential for %d planets:\n", n_test, n);

  /* Bounding box */
  bmin[0] = -1.0;
  bmin[1] = -1.0;
  bmin[2] = -1.0;
  bmax[0] = 1.0;
  bmax[1] = 1.0;
  bmax[2] = 1.0;

  /* Create random state */
  st = (state *)malloc(sizeof(state));
  st->n = n;
  st->y = (real **)malloc(sizeof(real *) * 3);
  st->y[0] = (real *)malloc(sizeof(real) * n);
  st->y[1] = (real *)malloc(sizeof(real) * n);
  st->y[2] = (real *)malloc(sizeof(real) * n);
  for (i = 0; i < n; i++)
  {
    st->y[0][i] = (bmax[0] - bmin[0]) * rand() / RAND_MAX + bmin[0];
    st->y[1][i] = (bmax[1] - bmin[1]) * rand() / RAND_MAX + bmin[1];
    st->y[2][i] = (bmax[2] - bmin[2]) * rand() / RAND_MAX + bmin[2];
  }

  /* Create Chebyshev interpolation points */
  in = new_chebyshev_interpolation(m);

  /* Setup cluster tree */
  s = setup_cluster(st, bmin, bmax, in, 16, 64);
  printf("Tree setup:\n");

  /* Check evaluation */
  x1 = (real *)malloc(sizeof(real) * n_test);
  x2 = (real *)malloc(sizeof(real) * n_test);
  x3 = (real *)malloc(sizeof(real) * n_test);
  phi1 = (real *)malloc(sizeof(real) * n_test);
  phi2 = (real *)malloc(sizeof(real) * n_test);

  for (i = 0; i < n_test; i++)
  {
    x1[i] = (bmax[0] - bmin[0]) * rand() / RAND_MAX + bmin[0];
    x2[i] = (bmax[1] - bmin[1]) * rand() / RAND_MAX + bmin[1];
    x3[i] = (bmax[2] - bmin[2]) * rand() / RAND_MAX + bmin[2];
  }

  for (i = 0; i < n_test; i++)
    phi1[i] = eval_cluster(x1[i], x2[i], x3[i], s, st, eta);
  printf("Tree evaluation:\n");

  for (i = 0; i < n_test; i++)
    phi2[i] = eval_potential(x1[i], x2[i], x3[i], st);
  printf("Direct evaluation:\n");

  maxerror = 0.0;
  for (i = 0; i < n_test; i++)
  {
    error = ABS(phi1[i] - phi2[i]) / ABS(phi2[i]);
    maxerror = (error > maxerror ? error : maxerror);
  }
  printf("Maximal relative error %.4e\n", maxerror);

  return 0;
}
