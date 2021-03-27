function S=chol_P(A)

% Finding the Cholesky factorization of a symmetric, positive definite 
% matrix A

[m,n] = size(A);
if m~=n
    error('The matrix A is not square')
end

for i=1:n
    for j=1:i-1
        if A(i,j) ~=A(j,i)
                error('The matrix A is not symmetric')     
        end
    end
end

s = eig(A);

for i=1:n
    if s(i)<= 0
        error('A is not positive definite')
    end
end

S= zeros(n,n);

for i=1:n
    S(i,i) = sqrt(A(i,i) - norm(S(1:i-1,i),2)^2);
    for j=i+1:n
        temp = A(i,j) - S(1:i-1,i)'*S(1:i-1,j);
        S(i,j) = temp ./ S(i,i);
    end
end

S2 = chol(A);

S-S2