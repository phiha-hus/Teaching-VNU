#ifndef LINALG_HEADER
#define LINALG_HEADER

#include "basic.h"

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

#endif // LINALG_HEADER