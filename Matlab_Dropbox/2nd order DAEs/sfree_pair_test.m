clear all; close all; clc

m = 3; n = 5;
m1 = 3; r = 2; r1 = 2;

P = rand(n,n);

M = rand(m,m) * [rand(r,n) ; zeros(m-r,n)] * P
N = rand(m1,m1) * [rand(r1,n) ; zeros(m1-r1,n)] * P

sfree_pair(M,N)