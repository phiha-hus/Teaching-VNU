#include "cluster.h"

/* ------------------------------------------------------------
 * Constructor for Chebyshev interpolation
 * ------------------------------------------------------------ */

interpolation *
new_chebyshev_interpolation(int m)
{
    interpolation *in;
    int i;

    in = (interpolation *)malloc(sizeof(interpolation));
    in->m = m;
    in->xi = (real *)malloc(sizeof(real) * (m + 1));
    for (i = 0; i <= m; i++)
        in->xi[i] = cos(M_PI * (m - i + 0.5) / (m + 1.0));

    return in;
}

real eval_lagrange(int m, const real *xi, int nu, real x)
{
    real num, den;
    int i;

    num = 1.0;
    den = 1.0;
    for (i = 0; i < nu; i++)
    {
        num *= x - xi[i];
        den *= xi[nu] - xi[i];
    }
    for (i = nu + 1; i <= m; i++)
    {
        num *= x - xi[i];
        den *= xi[nu] - xi[i];
    }

    return num / den;
}

/* ------------------------------------------------------------
 * Constructor for a cluster
 * ------------------------------------------------------------ */

pcluster new_cluster(int size, int *idx, const real *bmin, const real *bmax)
{
    pcluster s;

    s = (pcluster)malloc(sizeof(cluster));

    s->size = size;
    s->idx = idx;

    s->sons = 0;
    s->son = 0;

    s->bmin = (real *)malloc(sizeof(real) * 3);
    s->bmin[0] = bmin[0];
    s->bmin[1] = bmin[1];
    s->bmin[2] = bmin[2];

    s->bmax = (real *)malloc(sizeof(real) * 3);
    s->bmax[0] = bmax[0];
    s->bmax[1] = bmax[1];
    s->bmax[2] = bmax[2];

    s->m = 0;
    s->xi = 0;

    return s;
}

/* ------------------------------------------------------------
 * Split a cluster
 * ------------------------------------------------------------ */

void splitsub(int size, int *idx, const real **y,
              real *bmin, real *bmax,
              int k,
              int stride, int *sizes, int **idxs, real *son_bmin, real *son_bmax)
{
    real c, d;
    int i, j, h;

    if (k == 0)
    {
        /* Store size of cluster */
        sizes[0] = size;

        /* Store index set */
        idxs[0] = idx;

        /* Store bounding box of cluster */
        son_bmin[0] = bmin[0];
        son_bmin[1] = bmin[1];
        son_bmin[2] = bmin[2];

        son_bmax[0] = bmax[0];
        son_bmax[1] = bmax[1];
        son_bmax[2] = bmax[2];

        /* Nothing more to do */
        return;
    }
    else
    {
        /* Decrease dimension */
        k--;

        /* Compute midpoint of current dimension's interval */
        c = (bmax[k] + bmin[k]) * 0.5;

        /* Split index set */
        i = 0;
        j = size - 1;
        while (i <= j)
        {
            /* Scan forwards as long as points belong in the "left" cluster */
            while (i <= j && y[k][idx[i]] < c)
                i++;
            /* Now idx[0],...,idx[i-1] belong in the "left" cluster */

            /* Scan backwards as long as points belong in the "right" cluster */
            while (i <= j && y[k][idx[j]] >= c)
                j--;
            /* Now idx[j+1],...,idx[size-1] belong in the "right" cluster */

            if (i < j)
            {
                /* Now we have y[k][idx[i]] >= c and y[k][idx[j]] < c, so we
	 * swap the indices */
                h = idx[i];
                idx[i] = idx[j];
                idx[j] = h;

                i++;
                j--;
            }
        }

        /* Process "left" cluster */
        d = bmax[k];
        bmax[k] = c;
        splitsub(i, idx, y, bmin, bmax, k, 2 * stride,
                 sizes, idxs, son_bmin, son_bmax);
        bmax[k] = d;

        /* Process "right" cluster */
        d = bmin[k];
        bmin[k] = c;
        splitsub(size - i, idx + i, y, bmin, bmax, k, 2 * stride,
                 sizes + stride, idxs + stride, son_bmin + 3 * stride, son_bmax + 3 * stride);
        bmin[k] = d;
    }
}

void split_cluster(pcluster s, const real **y)
{
    int *sizes;
    int **idxs;
    real *son_bmin, *son_bmax;
    int sons;
    int i, j;

    /* Initialize arrays for sons */
    sizes = (int *)malloc(sizeof(int) * 8);
    idxs = (int **)malloc(sizeof(int *) * 8);
    son_bmin = (real *)malloc(sizeof(real) * 8 * 3);
    son_bmax = (real *)malloc(sizeof(real) * 8 * 3);

    /* Split cluster */
    splitsub(s->size, s->idx, y,
             s->bmin, s->bmax,
             3,
             1, sizes, idxs, son_bmin, son_bmax);

    /* Count sons */
    sons = 0;
    for (i = 0; i < 8; i++)
        if (sizes[i] > 0)
            sons++;

    /* Create sons */
    s->son = (cluster **)malloc(sizeof(cluster *) * sons);
    i = 0;
    for (j = 0; j < 8; j++)
        if (sizes[j] > 0)
        {
            s->son[i] = new_cluster(sizes[j], idxs[j], son_bmin + 3 * j, son_bmax + 3 * j);
            i++;
        }
    assert(i == sons);
    s->sons = sons;

    free(son_bmax);
    free(son_bmin);
    free(idxs);
    free(sizes);
}

/* ------------------------------------------------------------
 * Create a cluster tree
 * ------------------------------------------------------------ */

void setup(pcluster s, const state *st, const interpolation *in,
           int resolution, int maxdepth)
{
    int m = in->m;
    const real *xi = in->xi;
    real middle, radius;
    real lag1, lag2, lag3;
    int nu1, nu2, nu3, mu1, mu2, mu3;
    int i, i1, i2, i3, j1, j2, j3;

    /* Set up interpolation points */
    assert(s->xi == 0);
    s->m = m;
    s->xi = (real **)malloc(sizeof(real *) * 3);
    s->xi[0] = (real *)malloc(sizeof(real) * (m + 1));
    s->xi[1] = (real *)malloc(sizeof(real) * (m + 1));
    s->xi[2] = (real *)malloc(sizeof(real) * (m + 1));

    middle = (s->bmax[0] + s->bmin[0]) * 0.5;
    radius = (s->bmax[0] - s->bmin[0]) * 0.5;
    for (i = 0; i <= m; i++)
        s->xi[0][i] = middle + xi[i] * radius;

    middle = (s->bmax[1] + s->bmin[1]) * 0.5;
    radius = (s->bmax[1] - s->bmin[1]) * 0.5;
    for (i = 0; i <= m; i++)
        s->xi[1][i] = middle + xi[i] * radius;

    middle = (s->bmax[2] + s->bmin[2]) * 0.5;
    radius = (s->bmax[2] - s->bmin[2]) * 0.5;
    for (i = 0; i <= m; i++)
        s->xi[2][i] = middle + xi[i] * radius;

    /* Set up interpolation coefficients */
    s->z = (real *)malloc(sizeof(real) * (m + 1) * (m + 1) * (m + 1));
    for (i = 0; i < (m + 1) * (m + 1) * (m + 1); i++)
        s->z[i] = 0.0;

    if (s->size > resolution && maxdepth > 0)
    {
        /* Split cluster */
        split_cluster(s, (const real **)st->y);

        /* Set up clusters */
        for (i = 0; i < s->sons; i++)
            setup(s->son[i], st, in, resolution, maxdepth - 1);

        /* Transfer coefficients by re-interpolation */
        for (i = 0; i < s->sons; i++)
        {

            // TODO: Call to build_transfer

            for (nu1 = 0; nu1 <= m; nu1++)
            {
                j1 = nu1;
                for (mu1 = 0; mu1 <= m; mu1++)
                {
                    i1 = mu1;
                    lag1 = eval_lagrange(m, s->xi[0], nu1, s->son[i]->xi[0][mu1]);
                    for (nu2 = 0; nu2 <= m; nu2++)
                    {
                        j2 = j1 * (m + 1) + nu2;
                        for (mu2 = 0; mu2 <= m; mu2++)
                        {
                            i2 = i1 * (m + 1) + mu2;
                            lag2 = lag1 * eval_lagrange(m, s->xi[1], nu2, s->son[i]->xi[1][mu2]);
                            for (nu3 = 0; nu3 <= m; nu3++)
                            {
                                j3 = j2 * (m + 1) + nu3;
                                for (mu3 = 0; mu3 <= m; mu3++)
                                {
                                    i3 = i2 * (m + 1) + mu3;
                                    lag3 = lag2 * eval_lagrange(m, s->xi[2], nu3, s->son[i]->xi[2][mu3]);
                                    s->z[j3] += lag3 * s->son[i]->z[i3];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {

        // TODO: call to build_leaf_clusterbasis

        /* Evaluate Lagrange polynomials in planets' positions */
        for (i = 0; i < s->size; i++)
        {
            i1 = s->idx[i];
            for (nu1 = 0; nu1 <= m; nu1++)
            {
                j1 = nu1;
                lag1 = eval_lagrange(m, s->xi[0], nu1, st->y[0][i1]);
                for (nu2 = 0; nu2 <= m; nu2++)
                {
                    j2 = j1 * (m + 1) + nu2;
                    lag2 = lag1 * eval_lagrange(m, s->xi[1], nu2, st->y[1][i1]);
                    for (nu3 = 0; nu3 <= m; nu3++)
                    {
                        j3 = j2 * (m + 1) + nu3;
                        lag3 = lag2 * eval_lagrange(m, s->xi[2], nu3, st->y[2][i1]);

                        s->z[j3] += lag3;
                    }
                }
            }
        }
    }
}

pcluster setup_cluster(const state *st,
                       const double *bmin,
                       const double *bmax,
                       const interpolation *in,
                       int resolution, int maxdepth)
{
    int n = st->n;
    int *idx;
    int i;
    pcluster s;

    idx = (int *)malloc(sizeof(int) * n);
    for (i = 0; i < n; i++)
        idx[i] = i;

    s = new_cluster(n, idx, bmin, bmax);

    setup(s, st, in, resolution, maxdepth);

    return s;
}

real dist2_cluster(real x1, real x2, real x3, pcluster s)
{
    real dist2;

    dist2 = sqr(x1 < s->bmin[0] ? s->bmin[0] - x1 : s->bmax[0] < x1 ? x1 - s->bmax[0] : 0.0);
    dist2 += sqr(x2 < s->bmin[1] ? s->bmin[1] - x2 : s->bmax[1] < x2 ? x2 - s->bmax[1] : 0.0);
    dist2 += sqr(x3 < s->bmin[2] ? s->bmin[2] - x3 : s->bmax[2] < x3 ? x3 - s->bmax[2] : 0.0);

    return dist2;
}

real diam2_cluster(pcluster s)
{
    real diam2;

    diam2 = sqr(s->bmax[0] - s->bmin[0]);
    diam2 += sqr(s->bmax[1] - s->bmin[1]);
    diam2 += sqr(s->bmax[2] - s->bmin[2]);

    return diam2;
}