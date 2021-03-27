function [x,err,k]=penalty(fun,grad_fun,h,grad_h,...
g,grad_g,x0,tol,kmax,kmaxd,meth,varargin)
% PENALTY Constrained  optimization with penalty
%  [X,ERR,K]=PENALTY(FUN,GRAD_FUN,H,GRAD_H,...
%  G,GRAD_G,X0,TOL,KMAX,KMAXD,METH)
%  computes a local minimizer of the cost function
%  FUN under the constraints H=0 and G>=0, by the
%  penalty method. X0 is the initial point, TOL is
%  the tolerance for the stopping test, KMAX is the
%  maximum number of allowed iterations.
%  GRAD_FUN, GRAD_H, and GRAD_G contain the gradient
%  of FUN, H, and G, respectively. The variables
%  H, G, GRAD_H, and GRAD_G can be set to [], if they
%  are not present. The solution of the corresponding
%  unconstrained minimization problem is performed
%  by calling either Matlab FMINSEARCH function
%  (if METH=0) or DESCENT function (if METH>0).
%  When METH>0, KMAXD and METH contain respectively
%  the maximum number of allowed iterations for the
%  function DESCENT and the choice of the descent
%  directions. When METH>1
%  [X,ERR,K]=PENALTY(FUN,GRAD_FUN,H,GRAD_H,...
%  G,GRAD_G,X0,TOL,KMAX,KMAXD,METH, HESS)
%  is the correct calling instruction.
%  If METH=1 HESS is the function handle associated
%  with the Hessian is required, if METH=2 HESS is a
%  suitable approximation of the Hessian at the step 0.
xk=x0(:); alpha0=1;
if meth==1, hess=varargin{1};
elseif meth==2, hess=varargin{1};
else  hess=[]; end
if ~isempty(h), [nh,mh]=size(h(xk)); end
if ~isempty(g), [ng,mg]=size(g(xk)); else, ng=[]; end
err=tol+1; k=0;
alphak=alpha0; alphak2=alphak/2; told=.1;
while err>tol && k< kmax
P=@(x)Pf(x,fun,g,h,alphak2,ng);
grad_P=@(x)grad_Pf(x,grad_fun,h,g,...
                   grad_h,grad_g,alphak,ng);
if meth==0
options=optimset('TolX',told);
[x,err,kd]=fminsearch(P,xk,options);
err=norm(x-xk);
else
[x,err,kd]=descent(P,grad_P,xk,told,kmaxd,meth,hess);
err=norm(grad_P(x));
end
if kd<kmaxd, alphak=alphak*10; alphak2=alphak/2;
else alphak=alphak*1.5; alphak2=alphak/2; end
k=k+1; xk=x; told=max([tol,told/10]);
end
end % end of the function penalty
function y=Pf(x,fun,g,h,alphak2,ng)
y=fun(x);
if ~isempty(h), y=y+alphak2*sum((h(x)).^2); end
if ~isempty(g), G=g(x);
for j=1:ng, y=y+alphak2*max([-G(j),0])^2; end
end
end % end of function Pf
function y=grad_Pf(x,grad_fun,h,g,...
                  grad_h,grad_g,alphak,ng)
y=grad_fun(x);
if ~isempty(h), y=y+alphak*grad_h(x)*h(x); end
if ~isempty(g), G=g(x); Gg=grad_g(x);
for j=1:ng
if G(j)<0, y=y+alphak*Gg(:,j)*G(j); end
end, end
end % end of function grad_Pf
