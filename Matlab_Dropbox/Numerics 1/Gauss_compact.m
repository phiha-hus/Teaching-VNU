function Gauss_compact(A)

% This function provides the Gauss compact scheme of the matrix A

[~,n] = size(A);

C(:,1) = A(:,1);

D = zeros(n,n);

D(1,2:n) = A(1,2:n)/A(1,1);

for l = 2:n-1
  C(l:n,l) = A(l:n,l) - C(l:n,1:l-1)* D(1:l-1,l);
  D(l,l+1:n) = ( A(l,l+1:n)-C(l,1:l-1)*D(1:l-1,l+1:n) ) .* C(l:n-1,l:n-1);
end

C

D