clear all; close all; clc

A = [0 1 0 0; 0 0 -1 0; 0 0 0 1; 0 0 2 0]

B = [0 1 0 -1]'

C = [1 0 0 0]

D=[0]

pkg load signal

[N1,D1] = ss2tf(A,B,C,D,1)

[Z1,P1,K1] = tf2zp(N1,D1)

x0 = [1 1 1 1]

sys = ss(A,B,C,D)

################################################################################
# What Group 1 did not do
# @Copyright Phi Ha
################################################################################
# Plot the response  y for u = 0, x0 is given

x0 = [1 1 1 1]';
[Y,T,X] = initial(sys,x0)

figure(1); clf;
plot(T,X,T,Y) 
legend('x','y')
grid on


 