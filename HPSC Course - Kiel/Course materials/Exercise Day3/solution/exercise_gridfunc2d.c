
#include "basic.h"
#include "gridfunc2d.h"

int main(int argc, char const *argv[])
{
    pgridfunc2d u, b;
    real norm;
    uint n;

    printf(
        "\n"
        "##############################################################################\n"
        "#                            Exercise gridfunc2d                             #\n"
        "##############################################################################\n"
        "\n");

    // Number of gridpoints in each direction --> n^2 total inner points
    n = 100;

    u = new_gridfunc2d(n);
    b = new_gridfunc2d(n);

    clear_gridfunc2d(b);
    // Init the grid function 'u' with some reasonable data
    init_sine_gridfunc2d(u);
    write_gridfunc2d(u, "data/gridfunction_u.dat");

    // setup right-hand-side by computing b = Lu;
    addeval_5point_stencil_gridfunc2d(b, 1.0, u);

    // solve again for x
    solve_gridfunc2d(b);

    // write grid function to file
    write_gridfunc2d(b, "data/gridfunction_u_solve.dat");

    // Test if both grid functions are equal
    norm = diffnorm_gridfunc2d(u, b);
    printf("error: %.5e\n", norm);

    return 0;
}
