clear all; close all; clc
N = 1000;

f = @(n) 1/n;

for n = 1:N
    disp(strcat('u_',num2str(n),'= ',num2str(f(n))))    
end


%%
x = -3:0.01:3
y = -2*x.^2+5*x+1

plot(x,y,'-.','linewidth',2)
legend('-2*x.^2+5*x+1')
grid on

