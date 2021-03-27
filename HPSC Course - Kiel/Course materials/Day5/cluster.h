#ifndef CLUSTER_HEADER
#define CLUSTER_HEADER

#include "basic.h"

/* Representation of the state of the system */
typedef struct
{
    /* Number of planets */
    int n;

    /*  
    * Positions of planets
    * 
    * y[0] array of all x-coordinates for all masses
    * y[1] array of all y-coordinates for all masses
    * y[2] array of all z-coordinates for all masses
    * 
    * E.g. y[0][i] x-coordinate of mass i
    */
    real **y;
} state;

/* Representation of an interpolation scheme */
typedef struct
{
    /* Degree of interpolation */
    int m;

    /* Interpolation points in [-1,1] */
    real *xi;
} interpolation;

/* Representation of a cluster */
typedef struct _cluster cluster;
typedef cluster *pcluster;
struct _cluster
{
    /* Index set */
    int size;
    int *idx;

    /* Sons */
    int sons;
    cluster **son;

    /* Bounding box */
    real *bmin;
    real *bmax;

    /* Interpolation points */
    int m;
    real **xi;

    /* Interpolation coefficients */
    real *z;
};

/* ------------------------------------------------------------
 * Constructor for Chebyshev interpolation
 * ------------------------------------------------------------ */

/**
 * Create new chebyshev interpolation scheme of order 'm'
 */
interpolation *new_chebyshev_interpolation(int m);

/**
 * Evaluate Lagrange polynomial L_\nu of order 'm' at position
 * 'x' with interpolation points xi_0 , ... , x_m. 
 */
real eval_lagrange(int m, const real *xi, int nu, real x);

/* ------------------------------------------------------------
 * Constructor for a cluster
 * ------------------------------------------------------------ */

pcluster new_cluster(int size, int *idx, const real *bmin, const real *bmax);

/* ------------------------------------------------------------
 * Split a cluster
 * ------------------------------------------------------------ */

void splitsub(int size, int *idx, const real **y,
              real *bmin, real *bmax,
              int k,
              int stride, int *sizes, int **idxs, real *son_bmin, real *son_bmax);

void split_cluster(pcluster s, const real **y);

/* ------------------------------------------------------------
 * Create a cluster tree
 * ------------------------------------------------------------ */

void setup(cluster *s, const state *st, const interpolation *in,
           int resolution, int maxdepth);

pcluster setup_cluster(const state *st,
                       const double *bmin,
                       const double *bmax,
                       const interpolation *in,
                       int resolution, int maxdepth);

real dist2_cluster(real x1, real x2, real x3, pcluster s);

real diam2_cluster(pcluster s);

#endif // CLUSTER_HEADER
