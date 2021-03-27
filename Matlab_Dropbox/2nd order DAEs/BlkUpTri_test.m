clear all; close all; clc

m = 5
n = 5
rm = 2 
rn = 3 
rp = 4

P = rand(m,m); Q = rand(n,n);
M = P * [rand(rm,n); zeros(m-rm,n)] * Q;
N = P * [rand(rn,n); zeros(m-rn,n)] * Q;
P = P * [rand(rp,n); zeros(m-rp,n)] * Q;
tic
[r2,r1,r0,A2,A1,A0] = BlkUpTri(M,N,P);
toc