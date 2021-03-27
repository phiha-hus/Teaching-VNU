#ifndef ODE_HEADER
#define ODE_HEADER

#include "basic.h"
#include "linalg.h"

/**
 * Function datatype for the representation of the right-hand-side of 
 * ordinary differential equations:
 *    y'(t) = f(t, y(t))
 * 
 * y(t) is assumed to be a vector-valued function and is therefore stored
 * in 'yt'. 
 * 't' is current time parameter, the result of the function evaluation
 * will be stored within 'yt1'.
 * Additional data might be passed via 'data'.
 */
typedef void (*ode_func)(real t, pvector yt, pvector yt1, void *data);

/**
 * Perform one step of the explicit Euler method.
 * 
 * 't' is the current time parameter.
 * 'yt' is the current value of y(t).
 * 'f' the right-hand-side,
 * 'delta' timestep size.
 * 'yt1' Some auxiliary vector, might be used for the evaluation of 'f'.
 * 'data' Additional data needed by 'f'
 * 
 * After completing the function the new values y(t + delta) should be stored within
 * 'yt'.
 */
void euler_step(real t, pvector yt, ode_func f, real delta, pvector yt1, void *data);

/**
 * Perform one step of the Runge method.
 * 
 * 't' is the current time parameter.
 * 'yt' is the current value of y(t).
 * 'f' the right-hand-side,
 * 'delta' timestep size.
 * 'yt1' Some auxiliary vector, might be used for the evaluation of 'f'.
 * 'data' Additional data needed by 'f'
 * 
 * After completing the function the new values y(t + delta) should be stored within
 * 'yt'.
 */
void runge_step(real t, pvector yt, ode_func f, real delta, pvector yt1, void *data);

/**
 * Perform one step of the Leapfrog method.
 * 
 * 't' is the current time parameter.
 * 'yt' is the current value of y(t).
 * 'f' the right-hand-side,
 * 'delta' timestep size.
 * 'yt1' Some auxiliary vector, might be used for the evaluation of 'f'.
 * 'data' Additional data needed by 'f'
 * 
 * After completing the function the new values y(t + delta) should be stored within
 * 'yt'.
 */
void leapfrog_step(real t, pvector yt, ode_func f, real delta, pvector yt1, void *data);

/**
 * Perform one step of the Crank-Nicolson method.
 * 
 * 't' is the current time parameter.
 * 'yt' is the current value of y(t).
 * 'f' the right-hand-side,
 * 'delta' timestep size.
 * 'yt1' Some auxiliary vector, might be used for the evaluation of 'f'.
 * 'data' Additional data needed by 'f'
 * 
 * After completing the function the new values y(t + delta) should be stored within
 * 'yt'.
 */
void crank_nicolson_step(real t, pvector yt, ode_func f, real delta, pvector yt1, void *data);

#endif // ODE_HEADER