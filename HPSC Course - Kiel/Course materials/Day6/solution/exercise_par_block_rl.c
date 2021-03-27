
#include "basic.h"
#include "linalg.h"

/* ------------------------------------------------------------
 * Setup test matrices
 * ------------------------------------------------------------ */

pmatrix
new_diaghilbert_matrix(uint rows)
{

    pmatrix a;
    field *aa;
    field sum;
    uint lda;
    uint i, j;

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

/* ============================================================
 * Main program
 * ============================================================ */

int main(int argc, char **argv)
{

    pmatrix a, a_entry, a_block, a_par_block;
    pvector b, x_entry, x_block, x_par_block;
    pstopwatch sw = new_stopwatch();
    real err_abs, norm_b, t;
    uint n; // Number of rows and columns for the entire matrix
    uint m; // Number division in row and in column direction. NOT THE SIZE OF A BLOCK!!!
    uint i;

    printf("Testing matrix RL decomposition\n");

    for (n = 2048; n <= 2048; n *= 2)
    {

        for (m = 32; m <= 32; m *= 2)
        {
            printf("============================================================\n");

            // Setup test matrix
            printf("Setting up test matrix of size %d x %d:\n", n, n);
            start_stopwatch(sw);
            a = new_diaghilbert_matrix(n);
            t = stop_stopwatch(sw);
            printf("  %.3f ms\n", t * 1.0e3);
            printf("  %.3f KB\n", n * n * sizeof(double) / 1024.0);

            printf("============================================================\n");

            printf("Clone matrix:\n");
            start_stopwatch(sw);
            a_entry = clone_matrix(a);
            a_block = clone_matrix(a);
            a_par_block = clone_matrix(a);
            t = stop_stopwatch(sw);
            printf("  %.3f ms\n", t * 1.0e3);

            printf("============================================================\n");

            printf("Setting up test vector:\n");
            start_stopwatch(sw);
            b = new_vector(n);
            clear_vector(b);
            x_entry = new_vector(n);
            clear_vector(x_entry);
            x_block = new_vector(n);
            clear_vector(x_block);
            x_par_block = new_vector(n);
            clear_vector(x_par_block);

            // Setup test vector
            for (i = 0; i < n; i++)
            {
                x_entry->x[i] = x_block->x[i] = x_par_block->x[i] = 1.0 / (1.0 + i);
            }
            t = stop_stopwatch(sw);
            printf("  %.3f ms\n", t * 1.0e3);

            printf("============================================================\n");

            printf("Computing right-hand-side vector by standard MVM:\n");
            start_stopwatch(sw);
            // Compute matrix-vector-product by standard gemv
            gemv(false, n, n, 1.0, a->x, a->ld, x_entry->x, 1, b->x, 1);
            norm_b = nrm2(n, b->x, 1);
            t = stop_stopwatch(sw);
            printf("  %.3f ms\n", t * 1.0e3);

            printf("\n\n");

            printf("============================================================\n");

            printf("Compute RL-factorization, per-entry version:\n");
            start_stopwatch(sw);
            // Decompose A = R * L
            rldecomp(n, a_entry->x, n);
            t = stop_stopwatch(sw);
            printf("  %.3f ms\n", t * 1.0e3);

            printf("============================================================\n");

            printf("Evaluate factorization:\n");
            start_stopwatch(sw);
            // Evaluate x <-- (R * L) * x
            rleval(n, a_entry->x, n, x_entry->x, 1);
            t = stop_stopwatch(sw);
            printf("  %.3f ms\n", t * 1.0e3);

            printf("============================================================\n");

            printf("Compute error of factorization:\n");
            start_stopwatch(sw);
            // Compute error between b and the in-place product (R * L) * x
            axpy(n, -1.0, b->x, 1, x_entry->x, 1);
            err_abs = nrm2(n, x_entry->x, 1);
            printf("  Absolute max error: %.3e\n", err_abs);
            printf("  Relative max error: %.3e\n", err_abs / norm_b);
            t = stop_stopwatch(sw);
            printf("  %.3f ms\n", t * 1.0e3);

            printf("\n\n");

            printf("============================================================\n");

            printf("Compute RL-factorization, block version:\n");
            printf("Using %d x %d blocks --> blocksize = %.3f KB\n",
                   n / m, n / m, n / m * n / m * sizeof(double) / 1024.0);
            start_stopwatch(sw);
            // Decompose A = R * L
            block_rldecomp(n, a_block->x, n, m);
            t = stop_stopwatch(sw);
            printf("  %.3f ms\n", t * 1.0e3);

            printf("============================================================\n");

            printf("Evaluate factorization:\n");
            start_stopwatch(sw);
            // Evaluate x <-- (R * L) * x
            rleval(n, a_block->x, n, x_block->x, 1);
            t = stop_stopwatch(sw);
            printf("  %.3f ms\n", t * 1.0e3);

            printf("============================================================\n");

            printf("Compute error of factorization:\n");
            start_stopwatch(sw);
            // Compute error between b and the in-place product (R * L) * x
            axpy(n, -1.0, b->x, 1, x_block->x, 1);
            err_abs = nrm2(n, x_block->x, 1);
            printf("  Absolute max error: %.3e\n", err_abs);
            printf("  Relative max error: %.3e\n", err_abs / norm_b);
            t = stop_stopwatch(sw);
            printf("  %.3f ms\n", t * 1.0e3);

            printf("\n\n");

            printf("============================================================\n");

            printf("Compute RL-factorization, parallel block version:\n");
            printf("Using %d x %d blocks --> blocksize = %.3f KB\n",
                   n / m, n / m, n / m * n / m * sizeof(double) / 1024.0);
            start_stopwatch(sw);
            // Decompose A = R * L
            par_block_rldecomp(n, a_par_block->x, n, m);
            t = stop_stopwatch(sw);
            printf("  %.3f ms\n", t * 1.0e3);

            printf("============================================================\n");

            printf("Evaluate factorization:\n");
            start_stopwatch(sw);
            // Evaluate x <-- (R * L) * x
            rleval(n, a_par_block->x, n, x_par_block->x, 1);
            t = stop_stopwatch(sw);
            printf("  %.3f ms\n", t * 1.0e3);

            printf("============================================================\n");

            printf("Compute error of factorization:\n");
            start_stopwatch(sw);
            // Compute error between b and the in-place product (R * L) * x
            axpy(n, -1.0, b->x, 1, x_par_block->x, 1);
            err_abs = nrm2(n, x_par_block->x, 1);
            printf("  Absolute max error: %.3e\n", err_abs);
            printf("  Relative max error: %.3e\n", err_abs / norm_b);
            t = stop_stopwatch(sw);
            printf("  %.3f ms\n", t * 1.0e3);

            printf("============================================================\n");

            printf("\n\n");

            printf("Cleanup:\n");

            del_matrix(a);
            del_matrix(a_entry);
            del_matrix(a_block);
            del_matrix(a_par_block);
            del_vector(x_entry);
            del_vector(x_block);
            del_vector(x_par_block);
            del_vector(b);
        }
    }

    return EXIT_SUCCESS;
}