function [x,n] = newton(f,df,x0,nmax,tol,delta,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function newton is an implementation of Newton method for solving f(x)=0
% Form: [x,n] = newton(f,df,x0,nmax,tol,delta,varargin)
% f: function 
% df: the derivatice function
% x0 : starting point 
% nmax: maximal number of steps to be taken
% tol: maximal allowable error between two consecutive iterations
% delta: maximal allowable of f'(x0)
% x: the approximated root
% n: the number iterated steps
% Stopping Criteria: The method stops after nmax iterations or the difference 
% between two consecutive iterates is smaller than tol
% Copyright: Phi Ha, September 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use nargin only for advanced students 

disp("Please input both f(x) and f'(x)-could be in finite difference")

if nargin ==3
  nmax = 1000; tol = 1e-6; delta = 1e-12;
elseif (nargin ~= 3 | nargin ~= 6)
  error('Wrong number of input arguments')
end

if abs(df(x0)) < delta
  error("Derivative of function is nearly zero. Stops.")
end

x = x0; n = 0; diff = 1;

for i=1:nmax
  %if abs(y)<= tol       % Do not use this, since the number of significant ...
                         % ... figures in x is not guarantee
  if abs(diff) <= tol   
    disp('Newton method converges since |x_{n+1} - x_n| <= tol.')
    break
  else
    xnew = x - df(x)\f(x); % The division df(x)\f(x) is also applicable for matrix ... 
                          % ... so it is better than f(x)/df(x)
                          % eventhough they have the same meaning
    diff = x - xnew; 
    x = xnew;
    n = n + 1; 
  end
end

if (abs(diff) > tol  && n == nmax)
  disp("Newton method may diverges.")
end
