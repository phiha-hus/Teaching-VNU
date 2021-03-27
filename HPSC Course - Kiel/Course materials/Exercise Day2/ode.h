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

/**
 * Perform one step of the classic fourth-order Runge-Kutta method.
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
void rk_classic_step(real t, pvector yt, ode_func f, real delta, pvector yt1, void *data);

/**
 * 't' is the current time parameter.
 * 'yt' is the current value of y(t).
 * 'm' is the order of the method.
 * 'dy' is a pointer to an array of vectors of length m+1. each element is supposed to contain the prior evaluations of 'f'.
 *      This array is created and initialized by this method is return via the pointer.
 * 'w' A pointer to an array of weights used by this method. The array is created and initialized and returned by this function.
 * 'f' the right-hand-side,
 * 'delta' timestep size.
 * 'data' Additional data needed by 'f'
 */
void adams_bashforth_init(real t, pvector yt, uint m, pvector **dy, real **w, ode_func f, real delta, void *data);

/**
 * Perform one step of the Adams-Bashforth method of order 'm'.
 * 
 * 't' is the current time parameter.
 * 'yt' is the current value of y(t).
 * 'm' is the order of the method.
 * 'dy' is an array of vectors of length m+1. each element contains the prior evaluations of 'f'
 *      it is supposed to be overwritten cyclically.
 * 'w' An array of weights used by the method.
 * 'i' denotes the current position within the array 'dy'
 * 'f' the right-hand-side,
 * 'delta' timestep size.
 * 'yt1' Some auxiliary vector, might be used for the evaluation of 'f'.
 * 'data' Additional data needed by 'f'
 * 
 * After completing the function the new values y(t + delta) should be stored within
 * 'yt'.
 */
void adams_bashforth_step(real t, pvector yt, uint m, pvector *dy, real *w, uint i, ode_func f, real delta, void *data);

#endif // ODE_HEADER