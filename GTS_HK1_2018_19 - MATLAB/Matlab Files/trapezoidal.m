function [value,err,n]=trapezoidal(f,a,b,nmax,eps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function trapezoidal aims to find the integral over [a,b] of f 
%                                       using trapezoidal quadrature rule
% The idea is based on RECURSIVE TRAPEZOIDAL FORMULAR (Cheney-Kincaid 7ed.)
% Form : [value,err,n] = trapezoidal(f,a,b,nmax,ep)
% f: function 
% a,b: begin & end points
% nmax: maximal number of steps to be taken
% eps: maximal allowable error
% value: the approximated value of the integral
% err: absolute error of the integral
% n: the number iterated steps
% Stopping criteria: The iteration stops after maximal number of iteration
% or when the error err is smaller to the desired error eps
% Compare to the function trapz in MATLAV/OCTAVE
% Copyright: Phi Ha, November 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use nargin only for advanced students 
if nargin ==3
  nmax = 10; eps = 1e-4;
elseif nargin < 5
  error('Wrong number of input arguments')
end

% Main program starts here

% This is the R(0,0) value (Cheney/Kincaid, p.199)
h = b-a;
k = 1;                      % k = 2^n is the number of nodes on [a,b]
value = h .* (feval(f,a) + feval(f,b))./2 ; 

tic

for n = 1:nmax
 h = h/2;                   % Decrease the step by 2 equal to increase n by 2
 S = 0;
 k = 2 * k;                 % Double the amount of nodes without compute 2^n
 
for i = 1:2:(k-1)
 %for i = 1:2:(2^n-1)
   S = S + feval(f,a + h .* i);
 end
 
err =  - value/2 + h .* S;

if abs(err) <= eps
   %disp('Trapezoidal method converges with ')
   %n
   value = value + err;
   break
 else 
   value = value + err;
 end 
end

toc


