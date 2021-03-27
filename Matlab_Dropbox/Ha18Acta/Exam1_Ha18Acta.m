% Section 1: Example 1 in the Article Ha18Acta
% Section 2: Example: Smith form may dramatically increase the index and smoothness
% requirements to the consider system
% Section 4: Example 2 in the Article Ha18Acta, non-uniquely solvable
% systems

clear all; close all; clc

%%
% Section 1
syms x y

P = [0 -1 x;0 0 -1;0 1 0]
Q = [0 0 0;1 0 0;0 0 0]
% R = P - y*Q
% det(R)
% [U1,P1,V1] = smithForm(P)

U = [0 0 1;0 -1 0; 1 x 1]; V = [0 0 1;1 0 0;0 1 0];
Pnew = U*P*V
Qnew = U*Q*V

%%
% Section 2
syms x
P = [x.^2 x-1;x x.^2]
[U,V,Si] = smithForm(P)

Q = [0 x-1 x 0; x 0 0 x; x 0 1 0; 0 x 0 -1]
[U2,V2,Si2] = smithForm(Q)

%% 
% Section 3
clear all; clc;

syms x
P = [-1 1 x;1 0 x; x 0 -1]
[U,V,Si] = smithForm(P,x)

%U * P * V - Si
%simplify(ans)

Q = [0 0 0; 0 1 0; 1 0 0]
U*Q*V
simplify(ans)

%% 
% Section 4

clear all; clc;

syms x
P = [0 -1 0; 0 x -1; 0 0 0]
Q = [0 -1 0; 0 0 0; 0 0 1]

U1 = [-1 0 0; -x -1 0; 0 0 1]
V1 = [0 1 0;0 0 1;1 0 0]'
P1 = U1*P*V1
Q1 = U1*Q*V1

U2 = eye(3)
V2 = [0 0 1; 1 0 0; 0 1 0]'
P2 = U2*P1*V2
Q2 = U2*Q1*V2

V1*V2

V1'-V2
V1*V2
V2*V1



