function U = Poisson(f,g,n,tol,max_it)
% 
% function U = Poisson(f,g,n,tol,max_it)
%
% The five point scheme for solving the Dirichlet BVP of the
% Poisson equation on the unit square.
% Input
%      f: the right hand side function
%      g: the Dirichlet boundary value function
%      n: the number of sub-intervals of [0,1]
%      tol: relative error tolerance of the iterative solution;
%           default value: 10^(-5)
%      max_it: maximal number of iterations allowed;
%           default value: 10,000
% Output
%      U: the solution u_{ij}, i,j=1,...,n+1
%
% To use the program, the user must supply two m-files, say
% 'f.m' and 'g.m', to define the right hand side function f
% of the differential equation and the boundary value function g.
% The user should also choose a positive integer n for the 
% number of subintervals of [0,1].
% A sample call would be
%      U = Poisson('f','g',n,1e-8,1000)
% with 10^(-8) as the tolerance for the relative errors of the 
% iterative solution, and a maximal number of 1,000 iterations 
% is allowed for solving the finite difference system.
% It is also possible to use
%     U = Poisson('f','g',n,1e-8)
% then the maximal number of iterations is the default valus 10,000.
% If the default values 10^(-5) and 10,000 are to be used for the
% iteration relative error tolerance and maximal number of iterations,
% then one can simply use
%     U = Poisson('f','g',n)

% if max_it, or max_it and tol, are not provided, use the default values
if nargin < 5
	max_it=10000
end
if nargin < 4
	tol=1e-5
end
% compute some parameters
n1 = n+1;
h = 1/n;
h2 = h*h/4;
Fr = zeros(n,n);
Fr(2:n,2:n) = h2*feval(f,(1:n-1)*h,(1:n-1)*h);
% initialization
U = zeros(n1,n1);
% specify boundary conditions
U(1,1:n1) = feval(g,0,(0:n)*h);
U(n1,1:n1) = feval(g,1,(0:n)*h);
U(1:n1,1) = feval(g,(0:n)*h,0);
U(1:n1,n1) = feval(g,(0:n)*h,1);
tol=h2*tol;
% iteration
rel_err = 1;
itnum = 0;
while ((rel_err>tol) & (itnum<=max_it))
err = 0;
umax = 0;
  for j = 2:n
    for i = 2:n
      temp = (U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1))/4-Fr(i,j);
      dif = abs(temp-U(i,j));
      if (err <= dif)
        err = dif;
      end
      U(i,j) = temp;
      temp = abs(temp);
      if(umax <= temp)
        umax = temp;
      end
    end
  end
itnum = itnum+1;
rel_err = err/umax;
end
% plot the numerical solution
X = (0:h:n*h)';
Y = X;
surf(X,Y,U')
xlabel('x-axis')
ylabel('y-axis')
zlabel('The numerical solution')
s1 = sprintf('h=%6.4f',h)
title(['Solution: ',s1])

