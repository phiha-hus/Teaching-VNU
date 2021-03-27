% Example 5 of the Article Ha18_Acta

clear all; close all; clc

%%
E = [1 0 0; 0 0 1; 0 0 0];
A = [-2 0 0; 0 1 0; 0 0 1];
B = [1 0 0; 0 1/2 0; 0 0 1/2];

U = [1 2 2; 2 1 2; 2 2 1]
V = inv(U)

E = 10 * V*E*U; 
A = 10 * V*A*U;
B = 10 * V*B*U;

%%

E = [2 -4 -8; -8 -4 2;12 16 12]
A = [28 36 36; -12 -14 -24; -12 -24 -14]
B = [2 -6 -6; 2 9 4;2 4 9]

U * E * V
U * A * V
U * B * V
