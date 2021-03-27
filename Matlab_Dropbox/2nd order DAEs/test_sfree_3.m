clear all; close all; clc
% Test the solver sfree_triple for discrete-time multibody system
% Obtained by discretizing the (continuous time) 2nd order DAE
% Example 1. This example is taken from
% Example 2. This example is taken from

%% Example 1  
% Parameters of the continuous systems
m1=3; m2=2; m3=1; b1=5; b2=10; b3=15; b4=20; 
k1=5; k2=15; k3=15; k4=20;

% Matrix coefficients of the continuous time
M1 = diag([m1 m2 m3 0]);
D1 = [b1+b2 -b2 0 0; -b2 b2+b3 0 0; 0 -b3 b3+b4 -b4; 0 0 -b4 b4];
B1 = [k1+k2 -k2 0 0; -k2 k2+k3 0 0; 0 -k3 k3+k4 -k4; 0 0 -k4 k4];

% Matrix coefficients of the discrete time

h = 0.01;
A = M1;
B = -2*M1 + h*D1;
C = M1 - h*D1+ h^2*B1;

[ka,M,N,P]=sfree_triple(A,B,C)
%[ka,M,N,P]=sfree_triple(M1,D1,B1)

%% Example 2: Example 3 from Ascher-Petzold'93
nu = 1;
h=0.01; 

A = @(n)[1 0 0; 0 1 0; 0 0 0]
B = @(n)[-2+h*nu*nu 0 0; 0 -2+h*nu*nu 0; 0 0 0]
C = @(n)[1-h*nu*nu C1(n)*h; 0 1-h*nu*nu C1(n)*h; C1(n) 0]

[ka,M,N,P]=sfree_triple_LTV(0.5,A,B,C)




