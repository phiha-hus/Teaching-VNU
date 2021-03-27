
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

    /**
     * TODO: init some test data in 'u', 
     *       generate the right-hand-side 'b' with b = L*u,
     *       solve for u again.
     */

    // write grid function to file
    write_gridfunc2d(u, "data/gridfunction_u.dat");

    // write grid function to file
    write_gridfunc2d(b, "data/gridfunction_u_solve.dat");

    // Test if both grid functions are equal
    norm = diffnorm_gridfunc2d(u, b);
    printf("error: %.5e\n", norm);

    return 0;
}
