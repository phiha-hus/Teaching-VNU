%clear all; close all; clc
global h A Ad alpha

% System coefficients
h = 0.1; kp = 0.4451; K = 1.53; 
ki = 2.3046; tau = 0.0254;

A = - [ 0 1; 0 -1/tau]; 
Ad = - [0 0; -K*ki/tau -K*kp/tau] ; 

%[alpha] = find_alpha(h,A,Ad)
W0 = lambertw_matrix(0,-h*Ad*expm(h*A));
W = fsolve(@witer,W0);
S0 = W0/h - A;
lambda = eig(S0);
alpha = max(lambda)

[K1,alpha] = find_K1(h, K, kp, ki, tau)

K3 = find_K3(h, A, Ad,alpha)

format shorte, K3