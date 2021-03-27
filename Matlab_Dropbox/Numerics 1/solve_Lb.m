function x=solve_Lb(L,b)

% Matlab solver for lower triangular system
% U: upper triangular matrix
% column vector 

[~,n] = size(L);
x = zeros(n,1);

for i=1:n
    x(i) = (b(i) - L(1:i-1,i)' * x(1:i-1))/L(i,i);
end