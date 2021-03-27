function x=solve_Ub(U,b)

% Matlab solver for upper triangular system
% U: upper triangular matrix
% column vector 

[~,n] = size(U);
x = zeros(n,1);

for i=n:-1:1
    x(i) = (b(i) - U(i,i+1:n) * x(i+1:n))/U(i,i);
end