
#include "basic.h"
#include "gridfunc1d.h"

int main(int argc, char **argv)
{
  pgridfunc1d x, b, d, exact;
  pgridfunc1d *xgrid;
  pgridfunc1d *bgrid;
  pgridfunc1d *dgrid;

  real error;
  uint L, n;
  uint i;

  /**
   * Init the Grids
   */

  // Level of the finest grid

  printf("Creating grid hierarchy:\n");
  L = 5;

  xgrid = (pgridfunc1d *)allocmem((L + 1) * sizeof(pgridfunc1d));
  bgrid = (pgridfunc1d *)allocmem((L + 1) * sizeof(pgridfunc1d));
  dgrid = (pgridfunc1d *)allocmem((L + 1) * sizeof(pgridfunc1d));

  n = 1;
  for (i = 0; i <= L; i++)
  {
    n = 2 * n + 1;

    xgrid[i] = new_gridfunc1d(n);
    bgrid[i] = new_gridfunc1d(n);
    dgrid[i] = new_gridfunc1d(n);
    clear_gridfunc1d(xgrid[i]);
    clear_gridfunc1d(bgrid[i]);
    clear_gridfunc1d(dgrid[i]);
  }

  /**
   * Testing CG-Solver on finest grid
   */

  printf("Testing CG-Solver on finest grid:\n");

  x = new_gridfunc1d(n);
  b = new_gridfunc1d(n);
  d = new_gridfunc1d(n);
  exact = new_gridfunc1d(n);

  clear_gridfunc1d(x);
  clear_gridfunc1d(b);
  clear_gridfunc1d(d);

  init_sine_gridfunc1d(exact);

  clear_gridfunc1d(b);
  addeval_3point_stencil_gridfunc1d(b, 1.0, exact);
  cg_gridfunc1d(b);
  error = diffnorm_gridfunc1d(b, exact);
  printf("  cg error: %.3e\n", error);

  /**
   * Testing smoother (Jacobi or Gauss-Seidel) as solver
   */

  printf("Testing Smoother:\n");

  clear_gridfunc1d(b);
  addeval_3point_stencil_gridfunc1d(b, 1.0, exact);
  for (i = 0; i < 500; i++)
  {
    // rich_step_1d(0.5 / (n + 1) / (n + 1), b, x, d);
    // jac_step_1d(0.5, b, x, d);
    gs_step_1d(b, x);
    if (i % 20 == 0)
    {
      error = diffnorm_gridfunc1d(x, exact);
      printf("  Richardson / Jacobi / GS error:  %3d:  %.3e\n", i + 1, error);
    }
  }

  /**
   * Setup Multigrid solver
   */

  for (i = 0; i <= L; i++)
  {
    clear_gridfunc1d(xgrid[i]);
    clear_gridfunc1d(bgrid[i]);
    clear_gridfunc1d(dgrid[i]);
  }

  // Setup rhs on finest grid
  addeval_3point_stencil_gridfunc1d(bgrid[L], 1.0, exact);

  /**
   * Testing Multigrid
   */

  printf("Testing Multigrid solver:\n");
  for (i = 0; i < 20; i++)
  {
    mg_gridfunc1d(L, xgrid, bgrid, dgrid);
    error = diffnorm_gridfunc1d(xgrid[L], exact);
    printf("  multigrid error:  %3d:  %.3e\n", i + 1, error);
  }

  /**
   * cleanup
   */

  for (i = L; i-- > 0;)
  {
    del_gridfunc1d(dgrid[i]);
    del_gridfunc1d(bgrid[i]);
    del_gridfunc1d(xgrid[i]);
  }
  free(dgrid);
  free(bgrid);
  free(xgrid);

  del_gridfunc1d(x);
  del_gridfunc1d(b);
  del_gridfunc1d(d);

  return 0;
}
