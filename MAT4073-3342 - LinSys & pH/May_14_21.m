% Lecture on Controllability and Observability
clear all; close all; clc

%% Solving Sylvester equation AX + XB = C
% Finally - compare the result with built-in commands ...
... lyap and sylv in MATLAB
n = 3; 
A = rand(n,n) 
B = rand(n,n)
C = rand(n,n)

% Create a new system using Kronecker product kron and ...
...Vecoperator reshape
Anew = kron(eye(n),A) + kron(B',eye(n))
bnew = reshape(C,n^2,1)

% Solve the system
Xnew = Anew\bnew
X = reshape(Xnew,n,n)

disp('Test the solution - Wrong if error smaller than tol.')
tol = 1e-9 ;
if norm(A*X+X*B-C)>tol
    warning('Something wrong happens')
else
    disp('Everything is fine. The 2-norm of error is:') 
    norm(A*X+X*B-C)
end

disp('Compare with sylv in MATLAB R2015a')
norm(X - sylv(A,B,C))

%% Compute Gramians, Kalman canonical forms
n = 3; p = 1;  q = 2;

P = rand(n,n) ; 
A = P * diag( -rand(n,1) * 100 ) * inv(P) ; 

m = 2 ;  % m is the number of controllable modes
B = [rand(m,p) ; zeros(n-m,p)];

C = rand(q,n); D = rand(q,p);

sys = ss(A,B,C,D)

Wc = gram(sys,'c')

Wo = gram(sys,'o')


% ctrbf(sys) % Octave
[Abar,Bbar,Cbar,T,~] = ctrbf(A,B,C) % Matlab

%obsvf(sys)  % Octave
















