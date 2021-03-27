
#include "basic.h"
#include "linalg.h"

/* ------------------------------------------------------------
 * Setup test matrices
 * ------------------------------------------------------------ */

/* Simple 2 times 2 version */
pmatrix
new_2x2_matrix()
{

    pmatrix a;
    double *aa;

    a = new_matrix(2, 2);
    aa = a->x;

    aa[0] = -18.0;
    aa[1] = -15.0;
    aa[2] = 4.0;
    aa[3] = 3.0;

    return a;
}

/* Simple 4 times 4 version */
pmatrix
new_4x4_matrix()
{

    pmatrix a;
    double *aa;

    a = new_matrix(4, 4);
    aa = a->x;

    aa[0] = 2.0;
    aa[1] = 4.0;
    aa[2] = 6.0;
    aa[3] = -2.0;
    aa[4] = -1.0;
    aa[5] = 0.0;
    aa[6] = 1.0;
    aa[7] = -5.0;
    aa[8] = -3.0;
    aa[9] = -3.0;
    aa[10] = -1.0;
    aa[11] = 4.0;
    aa[12] = 3.0;
    aa[13] = 1.0;
    aa[14] = 6.0;
    aa[15] = 1.0;

    return a;
}

pmatrix
new_diaghilbert_matrix(int rows)
{

    pmatrix a;
    double *aa;
    double sum;
    int lda;
    int i, j;

    a = new_matrix(rows, rows);
    aa = a->x;
    lda = a->ld;

    for (j = 0; j < rows; j++)
    {
        sum = 1.0;
        for (i = 0; i < j; i++)
        {
            aa[i + j * lda] = 1.0 / (1.0 + i + j);
            sum += fabs(aa[i + j * lda]);
        }
        for (i = j + 1; i < rows; i++)
        {
            aa[i + j * lda] = 1.0 / (1.0 + i + j);
            sum += fabs(aa[i + j * lda]);
        }
        aa[j + j * lda] = sum;
    }

    return a;
}

pmatrix
new_hilbert_matrix(int rows)
{

    pmatrix a;
    double *aa;
    int lda;
    int i, j;

    a = new_matrix(rows, rows);
    aa = a->x;
    lda = a->ld;

    for (j = 0; j < rows; j++)
    {
        for (i = 0; i < rows; i++)
        {
            aa[i + j * lda] = 1.0 / (1.0 + i + j);
        }
    }

    return a;
}

pmatrix
build_hilbert_matrix()
{
    int n = 1024;

    return new_hilbert_matrix(n);
}

pmatrix
build_ddom_hilbert_matrix()
{
    int n = 1024;

    return new_diaghilbert_matrix(n);
}

int main(int argc, char const *argv[])
{
    pmatrix a;
    pvector b, x;
    double err_abs, norm_b;
    int n;
    int i;

    printf(
        "\n"
        "##############################################################################\n"
        "#                               Exercise BLAS                                #\n"
        "##############################################################################\n"
        "\n");

    // Setup test matrix
    //  a = new_2x2_matrix();
    //    a = new_4x4_matrix();
    //   a  = build_hilbert_matrix();
    a = build_ddom_hilbert_matrix();

    // Print the matrix values to stdout, for debugging purpose
    //    print_matrix(a);

    n = a->rows;

    b = new_vector(n);
    clear_vector(b);
    x = new_vector(n);
    clear_vector(x);

    // Setup test vector
    for (i = 0; i < n; i++)
    {
        x->x[i] = 1.0 / (1.0 + i);
    }

    // Compute matrix-vector-product by standard gemv
    gemv(false, n, n, 1.0, a->x, a->ld, x->x, 1, b->x, 1);

    // Decompose A = R * L
    decomp_ul(a);
    // Print the matrix values to stdout, for debugging purpose
    //    print_matrix(a);
    // Evaluate x <-- (R * L) * x
    eval_ul(a, x);

    // Compute error between b and the in-place product (R * L) * x
    norm_b = nrm2(n, b->x, 1);
    axpy(n, -1.0, b->x, 1, x->x, 1);
    err_abs = nrm2(n, x->x, 1);
    printf("  Absolute max error: %.3e\n", err_abs);
    printf("  Relative max error: %.3e\n", err_abs / norm_b);

    del_matrix(a);
    del_vector(x);
    del_vector(b);

    return 0;
}