clear all; close all; clc

N = 100;

n = 1:1:N;
x = 1./n;
y = sin(x);

figure(1); clf;
subplot(1,3,1)
plot(n,x,'r-.','linewidth',2)
title('Cac gia tri cua x_n = 1/n')
xlabel('n')
ylabel('x_n')
legend('x_n')
grid on

subplot(1,3,2)
plot(n,y,'r-.','linewidth',2)
title('Cac gia tri cua f(x_n) = sin(x_n)')
xlabel('n')
ylabel('f(x_n) = 1/sin(x_n)')
legend('f(x_n) = 1/sin(x_n)')
grid on

subplot(1,3,3)
plot(x,y,'r-.','linewidth',2)
title('Do thi cua ham f(x) = sin(x)')
xlabel('x_n')
ylabel('sin(x_n)')
legend('sin(x_n)')
grid on

