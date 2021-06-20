%%
% Magnitude scaling and Equivalence relation in OCTAVE
clear all; close all; clc

pkg load control 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vi du tuong tu voi BT3 cua bai giua ky 

% Magnitude scaling 

%help step   % See the syntax step
%-- Function File: [Y, T, X] = step (SYS)
%-- Function File: [Y, T, X] = step (SYS, T)
%-- Function File: [Y, T, X] = step (SYS, TFINAL)
%-- Function File: [Y, T, X] = step (SYS, TFINAL, DT)

A = [-0.1 2;0 -1]; B = [10;0.1]; C = [0.2 -1]; D = 0;
sys = ss(A,B,C,D) ;

disp('Find minimal realization of this system')
tol = 1e-9
sys=minreal(sys,tol)


figure(1); clf;
[y,t,x] = impulse(sys);
plot(t,x(:,1),t,x(:,2),t,y)
legend('x1','x2','y')
title('Plot the impulse response for the system')
grid on

figure(2); clf;
[y,t,x] = step(sys,10);
plot(t,x(:,1),t,x(:,2),t,y)
legend('x1','x2','y')
title('Plot the step response for the system')
grid on

M1 = max(abs(x(:,1)))  
M2 = max(abs(x(:,2)))
My = max(abs(y))

P = [My/M1 0 ; 0 My/M2] ; 
A = P * A * inv(P)
B = P * B
C = C * inv(P) 
sys = ss(A,B,C,D) ;

figure(3); clf;
[y,t,x] = step(sys,10);
plot(t,x(:,1),t,x(:,2),t,y)
legend('x1','x2','y')
title('Plot the step response for the system')
grid on

M1 = max(abs(x(:,1)))  
M2 = max(abs(x(:,2)))
My = max(abs(y))
disp('Max of an amplitude a for step input is: ')
10/My

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checking the equivalence and zero-equivalence of two systems
% Using Sylvester equation X = sylvester (A, B, C)

A1 = [2 1 2; 0 2 2; 0 0 1]
A2 = [2 1 1;0 2 1;0 0 -1]
P = sylvester(A1,-A2,zeros(3,3))

