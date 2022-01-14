clear all; close all; clc

A = [1 1;4 -20];
B = [1 -1;1 -1];
C = [1 0];
D = 0;

disp('Cac gia tri rieng cua A la:')
eig(A)

disp('Xay dung he su dung control system toolbox MATLAB: ')
sys = ss(A,B,C,D)

x0 = [1;1]; 
t0 = 0; tf = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y,t,x] = step(sys,10)

figure(1); clf;
semilogy(t,x(:,1),'r-',t,x(:,2),'b-.')
title('Trang thai x cua he')
xlabel('t')
ylabel('value')
legend('x1','x2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





