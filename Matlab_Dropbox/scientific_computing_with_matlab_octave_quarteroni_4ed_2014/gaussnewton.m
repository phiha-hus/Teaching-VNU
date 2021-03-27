function [x,err,iter]= gaussnewton(r,jr,x0,tol,...
    kmax,varargin)
%GAUSSNEWTON  Solves nonlinear least squares problems
%  [X,ERR,ITER]=GAUSSNEWTON(R,JR,X0,TOL,KMAX)
%  solves the nonlinear least squares by the Gauss-
%  Newton method. R and JR are the function handles
%  associated with the function R and its Jacobian,
%  respectively. X0 is the initial point for the se-
%  quence. TOL is the tolerance for the stopping test,
%  KMAX is the maximum number of allowed iterations.
err=tol+1; k=0; xk=x0(:);
rk=r(xk,varargin{:}); jrk=jr(xk,varargin{:});
while err>tol & k< kmax
[Q,R]=qr(jrk,0); dk=-R \ (Q'*rk);
xk1=xk+dk;
rk1=r(xk1,varargin{:});
jrk1=jr(xk1,varargin{:});
k=k+1;  err=norm(xk1-xk);
xk=xk1; rk=rk1; jrk=jrk1;
end
x=xk; iter=k;
if (k==kmax & err > tol)
 fprintf(['Gauss-Newton method stopped \n',...
 'without converging to the desired tolerance \n',...
 'because the maximum number of iterations was \n',...
 'reached\n']);
end
