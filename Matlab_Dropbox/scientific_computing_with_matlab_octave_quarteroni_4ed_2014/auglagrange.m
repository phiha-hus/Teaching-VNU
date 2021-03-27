function [x,err,k]=auglagrange(fun,grad_fun,h,grad_h,...
    x0,lambda0,tol,kmax,kmaxd,meth,varargin)
% AUGLAGRANGE Constrained optimization
%  [X,ERR,K]=AUGLAGRANGE(FUN,GRAD_FUN,H,GRAD_H,...
%  X0,LAMBDA0,TOL,KMAX,KMAXD,METH)
%  computes a local minimizer of the cost function
%  FUN under the constraints H=0, by the augmented
%  Lagrangian method. X0 is the initial point, TOL
%  the tolerance for the stopping test, KMAX the
%  maximum number of allowed iterations.
%  GRAD_FUN and GRAD_H contain the gradient of FUN
%  and H respectively. The solution of the associated
%  unconstrained minimization problem is performed
%  by calling either the Matlab FMINSEARCH function
%  (if METH=0) or the DESCENT function (if METH>0).
%  When METH>0, KMAXD and METH contain respectively
%  the maximum number of allowed iterations for the
%  function DESCENT and the choice of the descent
%  directions. When METH>1
%  [X,ERR,K]=AUGLAGRANGE(FUN,GRAD_FUN,H,GRAD_H,...
%  X0,LAMBDA0,TOL,KMAX,KMAXD,METHi, HESS)
%  is the correct calling instruction.
%  If METH=1 HESS is the function handle associated
%  with the Hessian is required, if METH=2 HESS is a
%  suitable approximation of the Hessian at the step 0.
alpha0=1;
if meth==1, hess=varargin{1};
elseif meth==2, hess=varargin{1};
else, hess=[]; end
err=tol+1; k=0; xk=x0(:); lambdak=lambda0(:);
if ~isempty(h), [nh,mh]=size(h(xk)); end
alphak=alpha0; alphak2=alphak/2; told=0.1;
while err>tol && k< kmax
L=@(x)Lf(x,fun,lambdak,alphak2,h);
grad_L=@(x)grad_Lf(x,grad_fun,lambdak,alphak,h,grad_h);
if meth==0
options=optimset('TolX',told);
[x,err,kd]=fminsearch(L,xk,options);
err=norm(x-xk);
else
[x,err,kd]=descent(L,grad_L,xk,told,kmaxd,meth,hess);
err=norm(grad_L(x));
end
lambdak=lambdak-alphak*h(x);
if kd<kmaxd, alphak=alphak*10; alphak2=alphak/2;
else alphak=alphak*1.5; alphak2=alphak/2; end
k=k+1; xk=x; told=max([tol,told/10]);
end
end % end auglagrange
function y=Lf(x,fun,lambdak,alphak2,h)
y=fun(x);
if ~isempty(h)
y=y-sum(lambdak'*h(x))+alphak2*sum((h(x)).^2); end
end % end function Lf
function y=grad_Lf(x,grad_fun,lambdak,alphak,h,grad_h)
y=grad_fun(x);
if ~isempty(h)
   y=y+grad_h(x)*(alphak*h(x)-lambdak); end
end % end function grad_Lf
