function [zero,res,niter]=newton(fun,dfun,x0,tol,...
                                 nmax,varargin)
%NEWTON Finds function zeros.
% ZERO=NEWTON(FUN,DFUN,X0,TOL,NMAX) tries to find the
% zero ZERO of the continuous and differentiable
% function FUN nearest to X0 using the Newton method.
% FUN and its derivative DFUN accept real scalar input
% x and return a real scalar value. If the search
% fails an error message is displayed. FUN and DFUN
% can be either inline functions or anonymous
% functions or they can be defined by external m-files.
% ZERO=NEWTON(FUN,DFUN,X0,TOL,NMAX,P1,P2,...) passes
% parameters P1,P2,... to functions: FUN(X,P1,P2,...)
% and DFUN(X,P1,P2,...).
% [ZERO,RES,NITER]=NEWTON(FUN,...) returns the value of
% the residual in ZERO and the iteration number at
% which ZERO was computed.
x = x0;
fx = feval(fun,x,varargin{:});
dfx = feval(dfun,x,varargin{:});
niter = 0; diff = tol+1;
while diff >= tol & niter < nmax
   niter = niter + 1;      diff = - fx/dfx;
   x = x + diff;           diff = abs(diff);
   fx = feval(fun,x,varargin{:});
   dfx = feval(dfun,x,varargin{:});
end
if (niter==nmax & diff > tol)
    fprintf(['Newton stopped without converging to',...
    ' the desired tolerance because the maximum\n ',...
    'number of iterations was reached\n']);
end
zero = x; res = fx;
return
