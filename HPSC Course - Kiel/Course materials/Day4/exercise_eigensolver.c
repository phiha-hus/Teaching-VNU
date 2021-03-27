

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
  /**
  * TODO:
  */
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
  /**
  * TODO:
  */
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
  /**
  * TODO:
  */
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

  /**
  * TODO:
  */

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
