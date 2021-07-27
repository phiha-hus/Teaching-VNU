function [K3] = find_K3(h, A, Ad,alpha)

h1 = 1e-2; t = 0 : h1 : h;

syms s

K3 = 0 ; 

for i = 1:length(t)
 J3 = int(norm( expm( (-A-alpha * eye(2))*t(i) + A*h) * Ad), s, 0, t(i));
 K3 = max(K3,J3);
end

K3


