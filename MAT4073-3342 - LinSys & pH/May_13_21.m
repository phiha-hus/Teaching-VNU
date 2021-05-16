% Realization problem
% function tf2ss
%
%% VD1  G(s)=(2s+1)/(s-1)
clear all; close all; clc
N = [2 1]; Q = [1 -1];
[A,B,C,D] = tf2ss(N,Q)  % Ket qua chinh xac la dang chinh tac dieu khien duoc

%% VD2 G(s) = [-(12s+6)/(3s+34)  (22s+23)/(3s+34)]
clear all; close all; clc
% Vi G co 2 cot nen phai lam tung cot mot
% Lap trinh thi 0 yeu cau he so dau cua Q = 1
% Chi can quy dong mau so tren moi cot, chu 0 can ca ma tran nhu ly thuyet
N1 = [-12 -6]; Q1 = [3 34];
[A1,B1,C1,D1] = tf2ss(N1,Q1);

N2 = [22 23]; Q2 = [3 34];
[A2,B2,C2,D2] = tf2ss(N2,Q2);

A = blkdiag(A1,A2) 
B = blkdiag(B1,B2)
C = [C1 C2]
D = [D1 D2]

%% VD4 cac em xem Example 7 (Chen)
% cung lam tren tung cot, quy dong mau so tren tung cot