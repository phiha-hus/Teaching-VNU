%% 
a=[-0.1 2;0 -1];b=[10;0.1];c=[0.2 -1];d=0;
[y,x,t]=step(a,b,c,d);
plot(t,y,t,x(:,1),t,x(:,2),'rd')
grid on
legend('y','x1','x2')

%%
% Practice find controllable canonical form
clear all; close all; clc

num1=[4 -2 -20;0 0 1];
den1=[2 5 2];         % First column of transfer function
[a1,b1,c1,d1] = tf2ss(num1,den1) ;

num2 = [0 3 6;0 1 1]  ; 
den2 = [1 4 4]; % Second column of transfer function
[a2,b2,c2,d2] = tf2ss(num2,den2) ;

A = blkdiag(a1,a2)
B = blkdiag(b1,b2)
C = [c1 c2]
D = [d1 d2]

%% 
% Practice controllability checking
% System 6.11 (inverted pendulum)
clc;

A = [0 1 0 0;0 0 -1 0;0 0 0 1; 0 0 5 0]
B = [0 1 0 -2]'

Kc = ctrb(A,B)
if rank(Kc)== size(A)
    disp('System is controllable')
else
    disp('Uncontrollable')
end

eig(A)

Wc = gram(A,B)   % Se o chay duoc vi he 0 on dinh








