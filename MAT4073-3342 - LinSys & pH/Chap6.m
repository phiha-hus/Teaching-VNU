%%
% Test some control functions in OCTAVE
clear all; close all; clc

n = 3; p = 1;  q = 2;

P = rand(n,n) ; 
A = P * diag( -rand(n,1) * 100 ) * inv(P) ; 

m = 2 ;  # m is the number of controllable modes
B = [rand(m,p) ; zeros(n-m,p)];

C = rand(q,n); D = rand(q,p);

sys = ss(A,B,C,D)

Wc = gram(sys,'c')

Wo = gram(sys,'o')


ctrbf(sys)

obsvf(sys)

