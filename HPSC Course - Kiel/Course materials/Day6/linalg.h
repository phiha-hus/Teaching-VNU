#ifndef LINALG_HEADER
#define LINALG_HEADER

#include "basic.h"
#include "miniblas.h"

struct _vector
{
    field *x;
    uint dim;
};

typedef struct _vector vector;
typedef vector *pvector;

pvector new_vector(uint dim);
void del_vector(pvector x);
void clear_vector(pvector x);
void add_vector(pvector x, field alpha, pvector y);

/****************************************************************************
 * dense matrices
 ***************************************************************************/

struct _matrix
{
    field *x;
    uint rows;
    uint cols;
    uint ld;
};

typedef struct _matrix matrix;
typedef matrix *pmatrix;

pmatrix new_matrix(uint rows, uint cols);
void del_matrix(pmatrix x);
void clear_matrix(pmatrix x);
pmatrix clone_matrix(pmatrix x);

/****************************************************************************
 * Band matrices
 ***************************************************************************/

struct _bandmatrix
{
    uint rows;
    uint cols;
    uint kl;
    uint ku;
    field *x;
};

typedef struct _bandmatrix bandmatrix;
typedef bandmatrix *pbandmatrix;

pbandmatrix new_bandmatrix(uint rows, uint cols, uint kl, uint ku);
void del_bandmatrix(pbandmatrix x);

void print_bandstorage_bandmatrix(pbandmatrix x);
void clear_bandmatrix(pbandmatrix x);

pbandmatrix setup_5point_stencil_bandmatrix(uint n);

void decomplr_bandmatrix(pbandmatrix x);

/****************************************************************************
 * linear algebra routines for dense matrices
 ***************************************************************************/

void leval(int n, field *a, int lda, field *x, int incx);
void reval(int n, field *a, int lda, field *x, int incx);
void rleval(int n, field *a, int lda, field *x, int incx);

void rldecomp(int n, field *a, int lda);
void block_rsolve(int n, int m, const field *r, int ldr,
                  field *x, int ldx);
void block_lsolve_trans(int n, int m, const field *l, int ldl,
                        field *x, int ldx);

void block_rldecomp(int n, field *a, int lda, uint blocks);
void par_block_rldecomp(int n, field *a, int lda, int blocks);

#endif // LINALG_HEADER