clear; clc;
% Exercise 3.8
disp('Exercise 3.8')
n = 21;
f = @(x) sin(2*pi.*x)
% x = [-1:1] - 21 node
x = -1 + 2/(n-1) * [0:n-1];
% create y to inter
y =f(x);

t = linspace(-1,1,100);

lagrange = polyfit(x,y,8);
lagrange = polyval(lagrange, t);

spi = spline(x,y,t);

ft = f(t);

plot(t, ft,'-r',t, lagrange,'+k',t, spi, '*b')
legend('f(x)', 'lagrange', 'spline')
grid on

% Ph?n v?i nhi?u ?ang b? sai!

%perturbedData = (-1).^(1 + [0:n-1])*10^-4 +f(x);
%agrangeWithPerturbed = polyval(polyfit(x, perturbedData,4), t);
%spiWithPerturbed = spline(x,perturbedData, t);
%plot(t, ft, t, lagrangeWithPerturbed,'+r', t, spiWithPerturbed,'*b')
%legend('f(x)', 'lagrange with perturbed', 'spline with perturbed')
%grid on


% Exercise 3.14
disp('Exercise 3.14')
flowRate = [0, 35, 0.125, 5, 0, 5, 1, 0.5, 0.125, 0];
fre = [0:10]*10
interpft(flowRate,100)
