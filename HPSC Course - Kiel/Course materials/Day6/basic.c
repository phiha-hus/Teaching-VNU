#include "basic.h"

struct _stopwatch
{
    double start;
    double current;
};

pstopwatch
new_stopwatch()
{
    pstopwatch sw;

    sw = (pstopwatch)allocmem(sizeof(stopwatch));

    return sw;
}

void del_stopwatch(pstopwatch sw)
{
    freemem(sw);
}

void start_stopwatch(pstopwatch sw)
{
    sw->start = omp_get_wtime();
}

real stop_stopwatch(pstopwatch sw)
{
    sw->current = omp_get_wtime();
    return (sw->current - sw->start);
}

real sqr(real x)
{
    return x * x;
};