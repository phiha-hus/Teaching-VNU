% Chapter 1. Solving scalar and system of nonlinear equations

clear all; close all; clc

%% matlab built-in function
% fzero doi it nhat 2 input, do la ham f va dieu kien ban dau x0. ...
% Dieu nay rat giong voi pp Newton.
% Chu y: Cac phuong phap nhu phan doi, lap don, day cung deu 0 su dung
% dieu kien ban dau x0

% VD1
f = @(x) x^2-1
x0 = 0.1;
fzero(f,x0)
fsolve(f,x0)

% VD2
f2 = @(x)[x^3-1; x^5-1]
%fzero(f2,x0)
fsolve(f2,x0)

% VD3

f3 = @(x)[x(1)^3-1; x(2)^3-8];

x0 = [0.3; 0.5];

fsolve(f3,x0)

%% BT3 trên l?p

f = @(x) x.^4 - 4 * x.^3 + 6 * x.^2 - 4 * x -16;

figure(); clf;
fplot(f,[-4 4])
grid on

x0 = -1;
fzero(f,x0)

x0 = 3;
fzero(f,x0)

