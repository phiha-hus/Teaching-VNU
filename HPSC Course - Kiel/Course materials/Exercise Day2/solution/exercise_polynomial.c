
#include "basic.h"

/**
 * Function type for general 1D-function 
 */
typedef real (*func)(real x);

/**
 * Compute the first derivative of a given function 'f' at 0, i.e. f'(0)
 * 
 * 'f' given function
 * 'h' step size of the approximation
 */
real df0(func f, real h)
{
    real val;

    assert(h > 0.0);

    val = (f(-2.0 * h) - 8.0 * f(-1.0 * h) + 8.0 * f(1.0 * h) - f(2.0 * h)) / (12.0 * h);

    return val;
}

real limit(real x)
{
    return (exp(x) - 1.0) / x;
}

/**
 * Compute the limit towards 0 of a given function 'f', i.e.  lim_{t->0} f(t)
 * 
 * 'f' given function
 * 'h' step size of the approximation
 */
real limg0(func f, real h)
{
    real val;

    assert(h > 0.0);

    val = 1.0 / 3.0 * f(h) - 2.0 * f(0.5 * h) + 8.0 / 3.0 * f(0.25 * h);

    return val;
}

int main(int argc, char **argv)
{

    func f = exp;
    func g = limit;
    real h;
    uint i;
    real val, error, old_error;

    printf("Evaluate first derivative f'(0):\n");

    h = 1.0;
    old_error = 1.0;
    for (i = 0; i < 10; i++)
    {
        val = df0(f, h);
        error = fabs(1.0 - val);
        printf("df(0) = %.15f   error: %.5e (x%.2f)\n", val, error, old_error / error);
        old_error = error;
        h *= 0.5;
    }

    printf("\n");

    printf("Evaluate lim_{t->0} g(t):\n");

    h = 1.0;
    old_error = 1.0;
    for (i = 0; i < 10; i++)
    {
        val = limg0(g, h);
        error = fabs(1.0 - val);
        printf("lim_{t->0} g(t) = %.15f   error: %.5e (x%.2f)\n", val, error, old_error / error);
        old_error = error;
        h *= 0.5;
    }

    return 0;
}