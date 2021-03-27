function [lambda,x,iter]=eigpower(A,tol,nmax,x0)
%EIGPOWER Computes the eigenvalue with maximum modulus
%  of a real matrix.
%  LAMBDA=EIGPOWER(A) computes with the power method
%  the eigenvalue of A of maximum modulus from an
%  initial guess which by default is an all one vector.
%  LAMBDA=EIGPOWER(A,TOL,NMAX,X0) uses an absolute
%  error tolerance TOL (the default is 1.e-6) and a
%  maximum number of iterations NMAX (the default is
%  100), starting from the initial vector X0.
%  [LAMBDA,V,ITER]=EIGPOWER(A,TOL,NMAX,X0) also returns
%  the eigenvector V such that A*V=LAMBDA*V and the
%  iteration number at which V was computed.
[n,m] = size(A);
if n ~= m, error('Only for square matrices'); end
if nargin == 1
   tol = 1.e-06;   x0 = ones(n,1);   nmax = 100;
end
x0 = x0/norm(x0);
pro = A*x0;
lambda = x0'*pro;
err = tol*abs(lambda) + 1;
iter = 0;
while err>tol*abs(lambda) & abs(lambda)~=0 & iter<=nmax
   x = pro;              x = x/norm(x);
   pro = A*x;            lambdanew = x'*pro;
   err = abs(lambdanew - lambda);
   lambda = lambdanew;   iter = iter + 1;
end
return
