% Chapter 5: Quarteroni - Linear systems - Giai gan dung he tuyen tinh Ax=b

%% Iterative methods: cac phuong phap lap

clear all; close all; clc

% A = sprandsym(n,density)    % Sparse symmetric random matrix

% Increase size of  m, n from 1e+3 upto 5e+3 and check processing time
m = 3e+3; n = m; density = 0.4;  

A = sprandn(m,n,density) ;     % Sparse normally distributed random matrix
b = rand(m,1);

issparse(A)

tic 
x = A\b;
toc

tic
[L,U] = lu(A);
y = L\b;
z = U\y;
toc

%% Test Jacobi and Gauss-Seidel

clear all; close all; clc
n = 1e+2;
A = rand(n,n); 
cond(A)
max(abs(eig(A)))

b = rand(n,1);
x0 = zeros(n,1); nmax = 1000; tol = 1e-12;

% Question: Why the solution is so strange?
[x,iter]= itermeth(A,b,x0,nmax,tol,'G');   % Test solver of Quarteroni
x'

x = A\b; 
x'

%% Compare \, PG, PCG
% Example 5.16

err1 = []; err2 = []; 
t1 = []; t2 = [];

for n=4:20
    An = hilb(n);
    cond(An);
    if cond(An)>1e+5
        warning('condition number is too big')
    end
    
    xn = ones(n,1);
    bn = An * xn;
    
    tic
    x = An\bn;
    err1 = [err1 norm(x-xn)/norm(xn)];
    t1 = [t1 toc];
    
    tic
    x = pcg(An,bn);
    err2 = [err2 norm(x-xn)/norm(xn)];
    t2 = [t2 toc];
    
end
n = 4:20
t1
t2
err1
err2
















