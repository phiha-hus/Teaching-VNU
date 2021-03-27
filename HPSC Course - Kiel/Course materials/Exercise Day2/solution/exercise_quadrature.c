#include "basic.h"

/**
 * function type for a 1D integrand f : R --> C.
 * 
 * 'x' is the evaluation point for the function.
 * 'data' is some additional data needed by the function.
 */
typedef field (*quad_func)(real x, void *data);

/**
 * Evaluate an integral by applying the Trapezoidal rule.
 * 
 * 'a' and 'b' are the integration limits.
 * 'f' is the function to be integrated.
 * 'f_data' additional data for the evaluation of 'f'
 */
field eval_trapezoidal_rule(real a, real b, quad_func f, void *f_data)
{
    return 0.5 * (b - a) * (f(a, f_data) + f(b, f_data));
}

/**
 * Evaluate an integral by applying the composite Simpson rule.
 * 
 * 'a' and 'b' are the integration limits.
 * 'intervals' is the number of subintervals. 'intervals' == 1 is equal to the Trapezoidal rule
 * 'f' is the function to be integrated.
 * 'f_data' additional data for the evaluation of 'f'
 */
field eval_composite_simpson_rule(real a, real b, uint intervals, quad_func f, void *f_data)
{
    uint i;
    real x, h;

    field val = 0.0;

    h = (b - a) / intervals;

    // initial left interval term
    x = a;
    val = f(x, f_data);
    for (i = 0; i < intervals - 1; i++)
    {
        //midpoint term
        x += 0.5 * h;
        val += 4.0 * f(x, f_data);

        //inner interval term
        x += 0.5 * h;
        val += 2.0 * f(x, f_data);
    }
    //last midpoint term
    x += 0.5 * h;
    val += 4.0 * f(x, f_data);

    //final right interval term
    x += 0.5 * h;
    val += f(x, f_data);

    val *= (b - a) / (6.0 * intervals);

    return val;
}

field normal_distribution(real x, void *data)
{
    real sigma = *((real *)data);

    return exp(-x * x / (2.0 * sigma * sigma));
}

int main(int argc, char const *argv[])
{

    quad_func f;
    real a, b;
    real sigma;
    real integralValue, realValue;
    real error, old_error;
    uint k;

    f = normal_distribution;
    a = -0.5;
    b = 0.5;
    sigma = 2.0;
    realValue = 9.896802673670128669414225441869e-01;

    printf("Exact integral value: %.30e\n", realValue);

    printf("Computing integral:\n");

    integralValue = eval_trapezoidal_rule(a, b, f, &sigma);
    printf("Approximate integral value by trapezoidal rule:\n"
           "\t%.15e (%.3e)\n",
           integralValue, fabs(integralValue - realValue) / realValue);

    printf("\n\n");

    old_error = 1.0;
    for (k = 0; k < 8; k++)
    {
        integralValue = eval_composite_simpson_rule(a, b, 1 << k, f, &sigma);
        error = fabs(integralValue - realValue) / realValue;
        printf("Approximate integral value by trapezoidal rule:\n"
               "\t%.15e (%.3e) (x%.2f)\n",
               integralValue, error, old_error / error);
        old_error = error;
    }

    return 0;
}
