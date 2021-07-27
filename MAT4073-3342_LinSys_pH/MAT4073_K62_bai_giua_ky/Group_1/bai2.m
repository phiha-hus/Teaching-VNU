clear all; close all; clc

pkg load signal

A = [0 1 0 0; 0 0 -1 0; 0 0 0 1; 0 0 2 0]
B = [0 1 0 -1]'
C = [1 0 0 0]
D=[0]

% State space 2 transfer function 
[N1,D1] = ss2tf(A,B,C,D,1)   # What can you conclude here?

% Find zeros and pole //   
% = 1 * (1*s^2-1)/(1*s^4-2s^2)
% Ly thuyet k = 1, zeros = +/- 1, poles = 0 , +- can(2)
[Z1,P1,K1] = tf2zp(N1,D1)    # What can you conclude here?

sys = ss(A,B,C,D)

################################################################################
# What Group 1 did not do
# @Copyright Phi Ha
################################################################################
# Plot the response  y for u = 0, x0 is given

x0 = [1 1 1 1]';

[Y,T,X] = initial(sys,x0) ;

figure(1); clf;
plot(T,Y,'r*') 
legend('y')
grid on
title('Phan hoi dau vao 0, x0 = [1 1 1 1]')


# Plot the response for impulse function and step function 
t_final = 10 
[Y1,T1,X1] = impulse(sys,t_final);
[Y2,T2,X2] = step(sys,t_final);

figure(2); clf;
plot(T1,Y1,'b-',T2,Y2,'r-.') 
legend('impulse','step')
grid on
title('Phan hoi dau vao xung (impulse) & buoc nhay (step)')




