% Magnitude scaling 

%help step   % See the syntax step
%-- Function File: [Y, T, X] = step (SYS)
%-- Function File: [Y, T, X] = step (SYS, T)
%-- Function File: [Y, T, X] = step (SYS, TFINAL)
%-- Function File: [Y, T, X] = step (SYS, TFINAL, DT)

A = [-2 0 0;1 0 1;0 -2 -2]; B = [1;0;1]; C = [1 -1 0]; D = 0;
sys = ss(A,B,C,D) ;

figure(2); clf;
[y,t,x] = step(sys,10);
plot(t,x(:,1),t,x(:,2),t,x(:,3),t,y)
legend('x1','x2','x3','y')
title('Plot the step response for the system')
grid on

M1 = max(abs(x(:,1)))  
M2 = max(abs(x(:,2)))
M3 = max(abs(x(:,3)))
My = max(abs(y))

P = [My/M1 0 0; 0 My/M2 0 ; 0 0 My/M3] ; 
A = P * A * inv(P)
B = P * B
C = C * inv(P) 
sys = ss(A,B,C,D) ;

figure(3); clf;
[y,t,x] = step(sys,10);
plot(t,x(:,1),t,x(:,2),t,x(:,3),t,y)
legend('x1','x2','x3','y')
title('Plot the step response for the system')
grid on

M1 = max(abs(x(:,1)))  
M2 = max(abs(x(:,2)))
M3 = max(abs(x(:,3)))
My = max(abs(y))
disp('Max of an amplitude a for step input is: ')
10/My