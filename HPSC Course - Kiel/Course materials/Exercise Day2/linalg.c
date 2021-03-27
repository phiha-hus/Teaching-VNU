#include "linalg.h"

pvector new_vector(uint dim)
{
    pvector x;

    x = allocmem(sizeof(vector));
    x->dim = dim;
    x->x = allocmem(dim * sizeof(field));

    return x;
}

void del_vector(pvector x)
{
    assert(x != NULL);

    freemem(x->x);
    freemem(x);
}

void clear_vector(pvector x)
{
    uint dim = x->dim;

    uint i;
    for (i = 0; i < dim; i++)
    {
        x->x[i] = 0.0;
    }
}

void add_vector(pvector x, field alpha, pvector y)
{
    uint dim = x->dim;

    uint i;

    assert(dim == y->dim);

    for (i = 0; i < dim; i++)
    {
        x->x[i] += alpha * y->x[i];
    }
}

pmatrix new_matrix(uint rows, uint cols)
{
    pmatrix x;

    x = allocmem(sizeof(matrix));
    x->rows = rows;
    x->cols = cols;
    x->ld = rows;
    x->x = allocmem(rows * cols * sizeof(field));

    return x;
}

void del_matrix(pmatrix x)
{

    assert(x != NULL);

    freemem(x->x);
    freemem(x);
}

void clear_matrix(pmatrix x)
{
    uint rows = x->rows;
    uint cols = x->cols;
    uint ld = x->ld;

    uint i, j;

    for (j = 0; j < cols; j++)
    {
        for (i = 0; i < rows; i++)
        {
            x->x[i + j * ld] = 0.0;
        }
    }
}