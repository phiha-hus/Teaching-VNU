function [E1,A2]=sfree_pair(E,A)

[m,n] = size(M);
[m1,n1] = size(A);
if (n~=n1) or (m~=m1)
    error('Sizes do not match')
end

Z2 = null(E');
%[U,S,V] = svd(E);
%Z2 = U(rank(E)+1:end,:)'
T2 = null(Z2'*A);
Y2 = orth(Z2'*A);
Z1 = orth(E*T2);

disp('The strangeness-free formulation is')
E1 = Z1'*E;
A2 = Y2' * Z2' * A;

if rank([E1;A2]) ~= n
disp('rank condition fails')
end