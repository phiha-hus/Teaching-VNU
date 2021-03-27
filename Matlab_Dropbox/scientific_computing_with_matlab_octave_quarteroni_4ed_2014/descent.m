function [x,err,iter]= descent(fun,grad,x0,tol,kmax,...
                               meth,varargin)
%DESCENT Descent method for optimization
%  [X,ERR,ITER]=DESCENT(FUN,GRAD,X0,TOL,KMAX,METH,HESS)
%  computes a local minimizer of function FUN by the
%  descent method with Newton directions (METH=1),
%  quasi-Newton directions (BFGS) (METH=2), gradient
%  directions (METH=3) or conjugate gradient directions
%  with Fletcher and Reeves beta_k (METH=41),
%  Polak and Ribiere beta_k (METH=42),
%  Hestenes and Stiefel beta_k (METH=43).
%  The steplength is computed by the backtracking
%  technique.  FUN, GRAD and HESS (the latter being
%  used only if METH=1) are function handles associated
%  with the cost function, its gradient and its Hessian
%  matrix, respectively. If METH=2, HESS is a matrix
%  approximating the Hessian of FUN at the initial
%  point X0. TOL is the tolerance for the stopping
%  test, while KMAX is the maximum allowed number of
%  iterations. The function backtrack is called inside.
if nargin>6
if meth==1, hess=varargin{1};
elseif meth==2, H=varargin{1}; end
end
err=tol+1; k=0; xk=x0(:); gk=grad(xk); dk=-gk;
eps2=sqrt(eps);
while err>tol & k< kmax
if meth==1;        H=hess(xk); dk=-H\gk; % Newton
elseif meth==2     dk=-H\gk;             % BFGS
elseif meth==3     dk=-gk;               % gradient
end
[xk1,alphak]= backtrack(fun,xk,gk,dk);
gk1=grad(xk1);
if meth==2 % BFGS update
  yk=gk1-gk; sk=xk1-xk; yks=yk'*sk;
  if yks> eps2*norm(sk)*norm(yk)
  Hs=H*sk;
  H=H+(yk*yk')/yks-(Hs*Hs')/(sk'*Hs);
  end
elseif meth>=40 % CG update
  if meth == 41
   betak=-(gk1'*gk1)/(gk'*gk); % FR
  elseif meth == 42
   betak=-(gk1'*(gk1-gk))/(gk'*gk); % PR
  elseif meth == 43
   betak=-(gk1'*(gk1-gk))/(dk'*(gk1-gk)); % HS
  end
  dk=-gk1-betak*dk;
end
xk=xk1; gk=gk1; k=k+1; xkt=xk1;
for i=1:length(xk1); xkt(i)=max([abs(xk1(i)),1]); end
err=norm((gk1.*xkt)/max([abs(fun(xk1)),1]),inf);
end
x=xk; iter=k;
if (k==kmax & err > tol)
 fprintf(['Descent method stopped \n',...
 'without converging to the desired tolerance \n',...
 'because the maximum number of iterations was \n',...
 'reached\n']);
end
