function [x,F,niter] = newtonsys(Ffun,Jfun,x0,tol,...
                                nmax, varargin)
%NEWTONSYS Finds a zero of a nonlinear system
% [ZERO,F,NITER]=NEWTONSYS(FFUN,JFUN,X0,TOL,NMAX)
% tries to find the vector ZERO, zero of a nonlinear
% system defined in FFUN with jacobian matrix defined
% in the function JFUN, nearest to the vector X0.
% The variable F returns the residual in ZERO
% while NITER returns the number of iterations needed
% to compute ZERO. FFUN and JFUN are MATLAB functions
% defined in M-files.
niter = 0; err = tol + 1; x = x0;
while err >= tol & niter < nmax
    J = feval(Jfun,x,varargin{:});
    F = feval(Ffun,x,varargin{:});
    delta = - J\F;
    x = x + delta;
    err = norm(delta);
    niter = niter + 1;
end
F = norm(feval(Ffun,x,varargin{:}));
if (niter==nmax & err> tol)
 fprintf([' Fails to converge within maximum ',...
         'number of iterations.\n',...
         'The iterate returned has relative ',...
         'residual %e\n'],F);
else
 fprintf(['The method converged at iteration ',...
         '%i with  residual %e\n'],niter,F);
end
return
