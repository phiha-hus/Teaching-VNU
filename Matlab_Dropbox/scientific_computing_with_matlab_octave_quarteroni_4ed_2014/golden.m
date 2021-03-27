function [xmin,fmin,iter]=golden(fun,a,b,tol,...
                                     kmax,varargin)
%GOLDEN Approximates a minimizer of 1D-functions
%  XMIN=GOLDEN(FUN,A,B,TOL,KMAX) approximates a
%  minimizer of the function FUN in the interval
%  [A,B] by the golden section method.
%  If the search fails, the program returns an
%  error message. FUN is either an anonymous
%  function, or an inline function, or a function
%  defined in a M-file.
%  XMIN=GOLDEN(FUN,A,B,TOL,KMAX,P1,P2,...) passes
%  parameters P1, P2,... to the function
%  FUN(X,P1,P2,...).
%  [XMIN,FMIN,ITER]= GOLDEN(FUN,...) returns
%  the value of FUN at XMIN and the number of
%  iterations required to compute XMIN.
phi=(1+sqrt(5))/2; phi1=1/phi; phi2=1/(phi+1);
c=a+phi2*(b-a); d=a+phi1*(b-a);
err=tol+1; k=0;
while err>tol & k< kmax
  if(fun(c) >= fun(d))
    a=c; c=d; d=a+phi1*(b-a);
  else
    b=d; d=c; c=a+phi2*(b-a);
  end
  k=k+1; err=abs(b-a)/(abs(c)+abs(d));
end
xmin=(a+b)/2; fmin=fun(xmin); iter=k;
if  (iter==kmax & err > tol)
 fprintf(['The golden section method stopped \n',...
 'without converging to the desired tolerance \n',...
 'because the maximum number of iterations was \n',...
 'reached\n']);
end
