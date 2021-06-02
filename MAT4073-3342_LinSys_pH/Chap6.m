%%
% Test some control functions in OCTAVE
clear all; close all; clc

pkg load control 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Teach class 4073
s = tf('s');
g = 1/(2*s^2+3*s+4);
h = c2d(g,0.1);   % Sampleing the system with h = 0.1

figure(1); clf;
step(h);
title ("Step response of a discretized PT2 transfer function");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 3; p = 1;  q = 1;
P = rand(n,n) ; 
A = P * diag( -rand(n,1) * 100 ) * inv(P) ; 
m = 2 ;  # m is the number of controllable modes
B = [rand(m,p) ; zeros(n-m,p)];
C = rand(q,n); D = rand(q,p);

sys = ss(A,B,C,D) ;

figure(2); clf;
[y,t,x] = step(sys,10);
plot(t,x(:,1),t,x(:,2),t,y)
legend('x1','x2','y')
title('Plot the step response for the system')
grid on



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Teach class MAT3342 
Wc = gram(sys,'c')
Wo = gram(sys,'ov')
ctrbf(sys)
obsvf(sys)

