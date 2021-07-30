function [x,n] = fix_iter_aitken(f,x0,nmax,eps,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function fix_iter_aitken is an implementation for the Fixed Point Iteration (phuong phap lap don)
% using Aitken extrapolation method
% f: function 
% x0 : starting point 
% nmax: maximal number of steps to be taken
% ep: maximal allowable error between two consecutive steps
% c: the approximated root
% n: the number iterated steps
% Stopping Criteria: The method stops after nmax iterations  
% or the difference between two consecutive iterates is smaller than eps
% Copyright: Phi Ha, September 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==2
   eps = 1e-4; nmax = 1e+3;
elseif nargin ==3
    disp('Specify both tolerance and maximum number of iterations')
elseif nargin > 4
    error('Number of input argumetns is wrong')
end

x = x0; y = f(x); diff = abs(y-x);
n = 0;

for i=1:nmax
    if diff <= eps
        disp('The Aitken extrapolation method converges.')
        break
    else
        z = f(y);
        xnew = (z .* x - y.^2)/(z - 2.* y + x); % Aitken extrapolation
        %xnew = f(x);   % Classical fixed_point_iteration
        diff = abs(xnew - x);
        n = n+1;
        x = xnew; 
        y = f(x);
    end
end

if (n == nmax && diff > eps)
    disp('Aitken extrapolation method diverges.')
end