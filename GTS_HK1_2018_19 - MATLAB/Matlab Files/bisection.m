function [c,err,n]=bisection(f,a,b,nmax,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function bisect aims to find the root of the function f using Bisection method
% Form : [c,err,n]=bisection(f,a,b,nmax,tol)
% f: function 
% a,b: begin & end points
% nmax: maximal number of steps to be taken
% tol: maximal allowable error
% c: the approximated root
% err: absolute error of the root
% n: the number iterated steps
% Stopping criteria: The iteration stops after maximal number of iteration
% or when the error err is smaller to the desired error tol
% Copyright: Phi Ha, September 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use nargin only for advanced students 
if nargin ==3
  nmax = 100; tol = 1e-4;
elseif nargin<5
  error('Wrong number of input arguments')
end

% Main program starts here

fa = f(a); fb = f(b);

if sign(fa) == sign(fb)
  error('Function has the same sign at a, b')
end

err = b-a;

for n=0:nmax
  err = err/2;
  c = a+err;
  fc = f(c);
  
  if err<tol
    disp('Bisection method converges.')        
    break
  end
  
  if sign(fa) ~= sign(fc)
    b=c;
    fb = fc;
  else 
    a=c;
    fa=fc;
  end
end
