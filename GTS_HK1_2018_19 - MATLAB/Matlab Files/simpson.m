function [value,err,n]=simpson(f,a,b,nmax,eps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function simpson aims to find the integral over [a,b] of f 
%                                       using Simpson rule
% The idea is based on ADAPTIVE SIMPSON QUADRATURE RULE (Cheney-Kincaid 7ed.)
% Form : [value,err,n] = simpson(f,a,b,nmax,ep)
% f: function 
% a,b: begin & end points
% nmax: maximal number of steps to be taken
% eps: maximal allowable error
% value: the approximated value of the integral
% err: absolute error of the integral
% n: the number iterated steps
% Stopping criteria: The iteration stops after maximal number of iteration
% or when the error err is smaller to the desired error eps
% Compare to the function quad in MATLAV/OCTAVE
% Copyright: Phi Ha, November 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use nargin only for advanced students 
if nargin ==3
  nmax = 10; eps = 1e-6;
elseif nargin<5
  error('Wrong number of input arguments')
end

% Main program starts here
f_end = feval(f,a) + feval(f,b) ; 
h = (b-a)/2;
k = 1;                      % k = 2^n is the number of nodes on [a,b]
I1 = feval(f,a+h); I2 = 0;
I = h/3 .* ( f_end + 4 * I1 ); 

for n = 1:nmax
 h = h/2;                   % Decrease the step by 2 equal to increase n by 2
 k = 2 * k;                 % Double the amount of nodes without compute 2^n
 
 I2 = I1 + I2;              % Update I2 new
 x = a+h:2*h:b-h;
 I1 = sum(f(x));
 
 I_new = h/3 * ( f_end + 4 * I1 + 2 * I2); 
 err = I_new - I;
 
  if abs(err) <= eps
     %disp('Trapezoidal method converges with ')
     %n
     value = I_new;
     break
  else 
     n = n+1;
  end

 value = I_new;
 
end