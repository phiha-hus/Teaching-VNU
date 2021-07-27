function [K3] = find_K3(h, A, Ad)

h=1;	A = 1;	Ad = 1;	
a = A; ad = Ad; 

W0 = lambertw_matrix(0,-h*Ad*expm(h*A));
alpha = real(W0/h - A);

K3 = norm(ad * (1 - expm(-a*h)) * expm(-alpha * h));



