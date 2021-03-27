% De ktra giua ky cho sinh vien 
% Bai 1: Ktra Gaussian elimination, LU decomposition, solving linear
% systems
% Bai 2: QR-decomposition 


%% De so 1 / Bai 1
clear all; close all; clc

A = [2 1 1 1;4 4 3 3;-2 1 0 4;8 8 10 7]
b = [6;18;4;41]
[L,U] = LU_dec(A)
x = A\b

[Q,R] = qr(A)
pause
clc

%% De so 1 / Bai 2
pause
clc



%% De so 2
clear all; close all; clc
 A = [3  1  2  4;-6 -1 -1 -8; -3 -3 -6 -2;-9 -5   -14   -12]
 b = [11;-17;-17;-45]
 [L,U] = LU_dec(A)
 x = A\b   
[Q,R] = qr(A)

