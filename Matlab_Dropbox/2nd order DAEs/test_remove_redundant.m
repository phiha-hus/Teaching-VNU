% Test set for verifying the removing redundancy Matlab function
% Example 1: Simplest one
% Example 2: is a very nice example since both Row1 and Row2 of P 
% do not belong to span(Row_1_of_Q,Row_2_of_Q) but the intersetion subspace
% of span(Row_1_of_P,Row_2_of_P) and span(Row_1_of_Q,Row_2_of_Q) 
% is of dimension 1.
% Example 3: Also like Example 2

% P = input('Input the matrix P = '); Q = input('Input the matrix Q = ')
% disp('Now we use the rows of Q to remove all hidden redundant equations in P')

%% Example 1

clear all; close all; clc
P = [1 0;0 0] 
Q = [0 0;1 0]
[S,Z1,Z2] = Remove_Redundant(P,Q)

%% Example 2 

clear all; close all; clc
P = [1 0 0;0 1 3;-1 1 3] 
Q = [1 2 3;0 0 1;0 0 1]
[S,Z1,Z2] = Remove_Redundant(P,Q)

%% Example 3

clear all; close all; clc
P = [1 1 2;0 3 5] 
Q = [2 0 0;3 1 3]
[S,Z1,Z2] = Remove_Redundant(P,Q)

%rank(S*Q)+rank(P)-rank([S*Q;P])
%Z1 * P + Z2 * Q
