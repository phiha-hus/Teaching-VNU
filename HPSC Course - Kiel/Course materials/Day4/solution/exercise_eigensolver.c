

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "basic.h"
#include "miniblas.h"
#include "linalg.h"

#define MAX_STEPS 100000

/*
 * Evaluates the 3-point-stencil matrix A given on the exerice sheet.
 * The parameter 'matrix' is not needed in this function, only the size of the
 * vectors 'x' and 'y' respectively is needed to compute the operation.
 * NO MEMORY ALLOCATION IS NEEDED!!!
 *
 * This function should be casted to type 'addeval_func'
 */
void addeval_1d_problem(double alpha, void *matrix, pvector x, pvector y)
{
  double *xx, *yx;

  int i, n;
  double h, diag, offdiag;

  (void)matrix;

  n = x->dim;
  xx = x->x;
  yx = y->x;

  assert(y->dim == n);
  assert(n >= 2);

  h = 1.0 / (n + 1);
  diag = 2.0 / (h * h);
  offdiag = -1.0 / (h * h);

  yx[0] += alpha * (diag * xx[0] + offdiag * xx[1]);
  for (i = 1; i < n - 1; i++)
  {
    yx[i] += alpha * (offdiag * xx[i - 1] + diag * xx[i] + offdiag * xx[i + 1]);
  }
  yx[n - 1] += alpha * (offdiag * xx[n - 2] + diag * xx[n - 1]);
}

/* ------------------------------------------------------------
 * Vector iteration
 * ------------------------------------------------------------ */

/*
 * Perform vector iteration in order to compute the biggest eigenvalue of a
 * matrix A. The matrix is given in terms of a callback function 'A' and some
 * additional parameters 'matrix'.
 * Within 'x' the eigenvector corresponding to the computed eigenvalue is stored.
 * 'eps' is the accuracy parameter. Use the right condition given by the lecture.
 * 'eigenvalue' The computed eigenvalue should be returned via this parameter.
 * 'resnorm' Also the residual norm should be returned via this parameter.
 *
 * Attention: this function should return the needed number of iterations.
 * It should be no more than 'MAX_STEPS' many steps.
 */
int power_iteration(addeval_func A, void *matrix, pvector x, double eps,
                    double *eigenvalue, double *resnorm)
{
  pvector y;
  double norm, rayleigh, residual;
  int i, n;

  rayleigh = 0.0;
  residual = 1.0;

  n = x->dim;
  y = new_vector(n);

  // start with random, unit vector 'x'
  random_vector(x);
  norm = nrm2(n, x->x, 1);
  scal(n, 1.0 / norm, x->x, 1);

  i = 0;
  while (i < MAX_STEPS && residual > eps * fabs(rayleigh))
  {
    scal(n, 0.0, y->x, 1); // set y <-- 0
    A(1.0, matrix, x, y);  // set y <-- Ax

    // Compute Rayleigh-quotient and residualnorm
    rayleigh = dot(n, x->x, 1, y->x, 1); // <x,x> == 1, because of scaling
    scal(n, -rayleigh, x->x, 1);
    axpy(n, 1.0, y->x, 1, x->x, 1);
    residual = nrm2(n, x->x, 1);

    copy_vector(y, x);
    norm = nrm2(n, x->x, 1); // 'norm' is current approximation of eigenvalue;

    scal(n, 1.0 / norm, x->x, 1);
    i++;
  }

  // write results back
  *eigenvalue = rayleigh;
  *resnorm = residual;

  // cleanup
  del_vector(y);

  return i; // Number of iterations needed.
}

/* ------------------------------------------------------------
 * inverse iteration
 * ------------------------------------------------------------ */

/*
 * Perform inverse iteration in order to compute the smallest eigenvalue of a
 * matrix A. The matrix is given in terms of a callback function 'A' and some
 * additional parameters 'matrix'.
 * Within 'x' the eigenvector corresponding to the computed eigenvalue is stored.
 * 'eps' is the accuracy parameter. Use the right condition given by the lecture.
 * 'eigenvalue' The computed eigenvalue should be returned via this parameter.
 * 'resnorm' Also the residual norm should be returned via this parameter.
 *
 * for solving a linear system Ax=b use the function 'cg_solve', which awaits
 * a callback-function for the matrix-vector-multiplication and a vector that
 * is used both as input 'b' and also as solution 'x'.
 *
 * Attention: this function should return the needed number of iterations.
 * It should be no more than 'MAX_STEPS' many steps.
 */
static int inverse_iteration(addeval_func A, void *matrix, pvector x,
                             double eps, double *eigenvalue, double *resnorm)
{
  pvector y;
  double norm, rayleigh, residual;
  int i, n;

  (void)matrix;

  rayleigh = 0.0;
  residual = 1.0;

  n = x->dim;
  y = new_vector(n);

  // start with random, unit vector 'x'
  random_vector(x);
  norm = nrm2(n, x->x, 1);
  scal(n, 1.0 / norm, x->x, 1);

  scal(n, 0.0, y->x, 1); // set y <-- 0
  A(1.0, matrix, x, y);  // set y <-- Ax

  // Compute Rayleigh-quotient and residualnorm
  rayleigh = dot(n, x->x, 1, y->x, 1); // <x,x> == 1, because of scaling
  axpy(n, -rayleigh, x->x, 1, y->x, 1);
  residual = nrm2(n, y->x, 1);

  i = 0;
  while (i < MAX_STEPS && residual > eps * fabs(rayleigh))
  {
    cg_solve(A, NULL, x);

    norm = nrm2(n, x->x, 1); // 'norm' is current approximation of eigenvalue;
    scal(n, 1.0 / norm, x->x, 1);

    scal(n, 0.0, y->x, 1); // set y <-- 0
    A(1.0, matrix, x, y);  // set y <-- Ax

    // Compute Rayleigh-quotient and residualnorm
    rayleigh = dot(n, x->x, 1, y->x, 1); // <x,x> == 1, because of scaling
    axpy(n, -rayleigh, x->x, 1, y->x, 1);
    residual = nrm2(n, y->x, 1);

    i++;
  }

  // write results back
  *eigenvalue = rayleigh;
  *resnorm = residual;

  // cleanup
  del_vector(y);

  return i; // Number of iterations needed.
}

/* ============================================================
 * Main program
 * ============================================================ */

int main(void)
{

  pvector x, v, w;
  double norm, lambda, eps;
  int n, iterations;

  n = 100;
  eps = 1.0e-10;

  srand(42);

  x = new_vector(n);
  random_vector(x);

  /* ------------------------------------------------------------
   * Vector iteration
   * ------------------------------------------------------------ */

  printf("Vector iteration\n");
  lambda = 0.0;
  norm = 0.0;
  iterations = 0;

  iterations = power_iteration((addeval_func)addeval_1d_problem, NULL, x,
                               eps, &lambda, &norm);

  printf("  Eigenvalue %.5e\n"
         "  Residual norm %.5e\n"
         "  Iterations %d\n",
         lambda, norm, iterations);

  /* ------------------------------------------------------------
   * Testing CG-Iteration
   * ------------------------------------------------------------ */

  printf("Testing CG-Iteration:\n");

  v = new_vector(n);
  w = new_vector(n);

  random_vector(v);
  clear_vector(w);
  addeval_1d_problem(1.0, NULL, v, w);
  norm = nrm2(n, v->x, 1);
  cg_solve(addeval_1d_problem, NULL, w);
  axpy(n, -1.0, v->x, 1, w->x, 1);
  printf("  rel. error: %.5e\n", nrm2(n, w->x, 1) / norm);

  del_vector(v);
  del_vector(w);

  /* ------------------------------------------------------------
   * Inverse iteration
   * ------------------------------------------------------------ */

  printf("Inverse iteration with cg\n");
  lambda = 0.0;
  norm = 0.0;
  iterations = 0;

  iterations = inverse_iteration((addeval_func)addeval_1d_problem, NULL, x,
                                 eps, &lambda, &norm);

  printf("  Eigenvalue %.5e\n"
         "  Residual norm %.5e\n"
         "  Iterations %d\n",
         lambda, norm, iterations);

  del_vector(x);

  return EXIT_SUCCESS;
}
