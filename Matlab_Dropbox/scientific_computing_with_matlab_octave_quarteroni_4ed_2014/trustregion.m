function [x,err,iter]= trustregion(fun,grad,x0,...
                delta0,tol,kmax,meth,hess)
%TRUSTREGION Trust region method for minimization
%  [X,ERR,ITER]=TRUSTREGION(FUN,GRAD,X0,TOL,KMAX,...
%  METH,HESS) computes a local minimizer of function
%  f by the trust region method. FUN and GRAD
%  (and HESS) are the function handles of the cost
%  function, its gradient (and its Hessian).
%  If METH=1, the Hessian HESS of f is used, otherwise
%  rank-one updates approximations of the Hessian are
%  built as in BFGS and the variable HESS is not requi-
%  red. X0 is the initial point for the sequence gene-
%  rated by the method. TOL is the tolerance for the
%  stopping test, KMAX is the maximum number of
%  iterations allowed.
delta=delta0; err=tol+1; k=0; mu=0.1;
eta1=0.25; eta2=0.75; gamma1=0.25; gamma2=2; deltam=5;
xk=x0(:);   gk=grad(xk);  eps2=sqrt(eps);
if meth==1 Hk=hess(xk); else Hk=eye(length(xk)); end
while err>tol & k< kmax
[s]=trustone(Hk,gk,delta);
rho=(fun(xk+s)-fun(xk))/(s'*gk+0.5*s'*Hk*s);
if rho> mu, xk1=xk+s; else, xk1=xk; end
if rho<eta1
     delta=gamma1*delta;
elseif rho> eta2 & abs(norm(s)-delta)<sqrt(eps)
    delta=min([gamma2*delta,deltam]);
end
gk1=grad(xk1);
err=norm((gk1.*xk1)/max([abs(fun(xk1)),1]),inf);
if meth==1 % Newton
   xk=xk1; gk=gk1; Hk=hess(xk);
else       % quasiNewton
  gk1=grad(xk1); yk=gk1-gk; sk=xk1-xk;
  yks=yk'*sk;
  if yks> eps2*norm(sk)*norm(yk)
  Hs=Hk*sk;
  Hk=Hk+(yk*yk')/yks-(Hs*Hs')/(sk'*Hs);
  end
  xk=xk1; gk=gk1;
end
k=k+1;
end
x=xk; iter=k;
if (k==kmax & err > tol)
 fprintf(['The trust region method stopped \n',...
 'without converging to the desired tolerance \n',...
 'because the maximum number of iterations was \n',...
 'reached\n']);
end
end

function [s]=trustone(Hk,gk,delta)
s=-Hk\gk; d = eigs(Hk,1,'sa');
if norm(s)>delta | d<0
lambda=abs(2*d); I=eye(size(Hk));
for l=1:3
R=chol(Hk+lambda*I);
s=-R \ (R'\gk); q=R'\s;
lambda=lambda+(s'*s)/(q'*q)*(norm(s)-delta)/delta;
if lambda< -d, lambda=abs(lambda*2); end
end
end
end
