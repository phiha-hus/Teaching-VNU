% Test functions for Example 1.3
% Notice that our equations are different from
% Duan's equations
clear all; close all; clc

h=1; A = -1; Ad = -1;

W0 = lambertw_matrix(0,h*Ad*expm(-h*A)) ;
W = fsolve(@witer,W0) ;
alpha = real(W0/h + A)


[K1] = find_K1(h, A, alpha) 

[K2] = find_K2(h,A,Ad,alpha)

[K3] = find_K3(h, A, Ad)

[K4] = find_K4(h,A,Ad,alpha)