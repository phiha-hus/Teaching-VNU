function [x,n] = fix_iter(f,x0,nmax,eps,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function fix_iter is an implementation for the Fixed Point Iteration (phuong phap lap don)
% Form : [x,n] = fix_iter(f,x0,nmax,eps)
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

% Use nargin only for advanced students 
if nargin ==2
  nmax = 100; eps = 1e-4;
elseif nargin<4
  error('Wrong number of input arguments')
end

x = x0; n = 0;

tic
for i = 1:nmax
  y = f(x); diff = abs(y-x);
  if diff <= eps    
    disp('Fixpoint iteration converges.') 
    break   
  else
     x = y;
     n = n+1;
  end 
end
toc

% tic
% for i=1:nmax
%     if abs(x-f(x))<= eps   % Wrong condition here, since x-f(x)
%         break
%     else
%         x = f(x);
%         n = n + 1;
%     end
% end
% toc

if (n == nmax && diff > eps)
    disp('Fix point iteration fails to converges.')
end

