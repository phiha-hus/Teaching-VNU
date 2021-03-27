function [x,alphak]= backtrack(fun,xk,gk,dk,varargin)
%BACKTRACK Backtracking strategy for line search.
%  [X,ALPHAK] = BACKTRACK(FUN,XK,GK,DK) computes the
%  new point x_{k+1}=x_k+alpha_k d_k, where alpha_k
%  is determined by the backtracking technique
%  with sigma=1.e-4 and rho=1/4.
%  [X,ALPHAK] = BACKTRACK(FUN,XK,GK,DK,SIGMA,RHO)
%  allows to specify the parameters sigma and rho.
%  Tipically 1.e-4<sigma<0.1 and 1/10< rho <1/2.
%  FUN is the function handle associated with the cost
%  function. XK, GK, and DK contain respectively the
%  point x_k, the gradient of f at x_k and the
%  descent direction d_k.
if nargin==4
    sigma=1.e-4; rho=1/4;
else
    sigma=varargin{1}; rho=varargin{2};
end
alphamin=1.e-5; % minimum value allowed for alpha_k
alphak = 1; fk = fun(xk);
k=0; x=xk+alphak*dk;
while fun(x)>fk+sigma*alphak*gk'*dk & alphak>alphamin
    alphak = alphak*rho;
    x = xk+alphak*dk; k = k+1;
end
