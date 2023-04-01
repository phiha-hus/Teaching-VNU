function [W info]=lambertwmatrix(k,A)

% W=LAMBERTWMATRIX(k,A) computes the determination k of the Lambert W function of A
%   k: branch of the Lambert W function
%   A: matrix whose Lambert W is required
%
% The algorithm is based on the Newton iteration with a heuristic stategy for
% the initial value based depending on the branch.
% The initial value is obtained truncating some series expansion of the 
% Lambert W function around the branch points and infinity.
% B. Iannazzo 2014
%
% Reference:
% M. Fasi, N. J. Higham and B. Iannazzo, An algorithm for the matrix 
% Lambert W function, MIMS Eprint 2014.58, November 2014

if size(k)~=[1 1]
  % control against common error
  error('Error: correct usage is lambertwmatrix(branch, MATRIX)')
end

% these parameter are obtained using a plot of the convergence regions
r0=Inf; % just for information
if (k==0)
  c0 = 1/2; % center of the disk for k=0
  r0 = 3/2; % default radius of the disk for k=0
  r0m = 1.35; r0M = 1.60; % smallest and largest radia of the disk for k=0
end
if (abs(k)==1)
  c0 = -1/2; % center of the disk for k=1,-1
  r0 = 1/3; % default radius of the disk for k=1,-1
  r0m = 0.25; r0M = 0.40; % smallest and largest radia of the disk for k=1,-1
end

n = size(A,1);

% first compute the Schur form of A to work on the triangular matrix Rt
[Qt Rt] = schur(A);
Q = Qt; 
R = Rt;

if (abs(k)<2)
  % find the largest gap to choose the splitting circle
  oeig = ordeig(Rt);
  woeig = oeig; % create a working copy of oeig 
  if (k==1)
    % leave alone eigenvalues with non negative imaginary part
    woeig(imag(woeig)>=0)=Inf;
  end
  if (k==-1)
    % leave alone eigenvalues with negative imaginary part
    woeig(imag(woeig)<0)=Inf;
  end
  [v s] = sort(abs(woeig-c0));
  v = [0;v;Inf]; % add c0 and Inf to simplify the procedure
  a = sum(v<r0m);
  b = sum(v>=r0m & v<=r0M); 
  c = sum(v>r0M);
  if (b~=0)
    % some eigenvalue in the annulus
    gap = [min(v(a+1)-v(a),2*(v(a+1)-r0m));...
           v(a+2:a+b)-v(a+1:a+b-1);...
           min(2*(r0M-v(a+b)),v(a+b+1)-v(a+b))];
    [m mind] = max(gap);  
    if (mind<=b)
      r0 = v(a+mind)-m/2;
    else
      r0 = r0M;
    end
  end
end

% we order the eigenvalues such that the resulting matrix is 2x2 block
% and the upper left block collect the eigenvalues for which the
% expansion 1 (around -1/e) is used, and the lower right block the ones
% for which the expansion 2 (around infty) is us
% the selection is made through a region (circle or half-circle) which is 
% different for k=0,1,-1
n1=0;
if (k == 0)
  select = abs(oeig-c0)<r0;
  [Q R] = ordschur(Qt,Rt,select);
  n1 = sum(select);
end
if (k == 1) 
  select = abs(oeig-c0)<r0 & imag(oeig)<0;
  [Q R] = ordschur(Qt,Rt,select);
  n1 = sum(select);
end
if (k == -1) 
  select = abs(oeig-c0)<r0 & imag(oeig)>=0;
  [Q R] = ordschur(Qt,Rt,select);
  n1 = sum(select);
end

% blocking, notice that R21 is 0
R11 = R(1:n1,1:n1);
R12 = R(1:n1,n1+1:n);
R22 = R(n1+1:n,n1+1:n);

% if the blocking is trivial use just one expansion
it1=0;it2=0;
if (n1 == 0)
  [W it2] = newton(R22,k,2);
end
if (n1 == n)
  [W it1] = newton(R11,k,1);  
end
% if the blocking is not trivial, compute the two diagonal blocks and then
% the upper right block
if (n1 > 0 && n1 < n)
  [W11 it1] = newton(R11,k,1);
  [W22 it2] = newton(R22,k,2);
  % block (1,2) of the Lambert
  C = R12*W22-W11*R12;
  W12 = sylvester(R11,-R22,-C);
  W = [W11 W12;zeros(n-n1,n1) W22];  
end

% revert the similarity which had put A in triangular form to get W(A)
W = Q*W*Q';

% return informations
info = sprintf('The branch is %d, number of eigs inside is %d\nNewton''s steps are %d and %d, respectively\nThe radius of splitting is %f',k,n1,it1,it2,r0);

end

function [W it]=newton(A,k,type)

% apply the Newton method to A with a starting  value depending on k and type
% k is the branch
% type is 1 for using asymptotic at -1/e
% type is 2 for using asymptotic at infinity

maxit=20;

n=size(A,1);
I=eye(n);

% choose the initial value
if (type == 1)
  % two terms of the expansion around -1/e
  X{1} = (1-2*abs(k))*sqrtm(2*exp(1)*A+2*eye(n)) - eye(n);
  H{1} = (A*expm(-X{1})-X{1})/(X{1}+I);
else 
  % three terms of the expansion around infty (ans 0 for |k|>0)
  L1 = logm(A) + 2*pi*i*k*I;
  L2 = logm(L1);
  L3 = L2/L1;
  X{1} = L1 - L2 + L3;
  H{1} = (L1*(expm(-L3) - I) + L2 - L3)/(X{1} + I);
end
X{2} = X{1} + H{1};

% Start the Newton iteration in its stable variant
it=2;
flag=0;
for k=3:maxit
  H{k-1} = ( ((X{k-2}+I)*H{k-2} + X{k-2} ) ...
               *expm(-H{k-2}) - X{k-1} ) / (X{k-1} + I);
  X{k} = X{k-1} + H{k-1};
  it = it + 1;
  if (flag==1) || norm(H{k-1})/norm(X{k}) < 1e-15  
    break;
  end
  if norm(H{k-1})/norm(X{k}) < 1e-8
    flag = 1;
  end
end

W = X{it};

end

function X = sylvester(A,B,Q)
% X=SYLVESTER(A,B,Q) solves the Sylvester equation AX + XB = Q 
% by using the Bartels and Stewart algorithm based on the complex
% Schur decomposition
%    A, B, Q: matrix coefficients
%    X : solution of AX + XB = Q 
[m,n] = size(Q);
[U,A1] = schur(A,'complex');
[V,B1] = schur(B.','complex');
Q1 = U'*Q*conj(V);
X = zeros(m,n);
X(:,n) = (A1 + B1(n,n)*eye(m))\Q1(:,n);
for i = n-1:-1:1
    v = Q1(:,i) - X(:,i+1:n)*B1(i,i+1:n).';
    X(:,i) = (A1 + B1(i,i)*eye(m))\v;
end
X = U*X*V.';

end
