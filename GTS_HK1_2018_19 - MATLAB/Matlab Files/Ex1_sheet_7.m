clear all; close all; clc

f = @(x) 1./(1+x.^2);
x = linspace(-5,5,21); y = f(x);

%v = div_diff(x,y)
z = linspace(-5,5,41); ynew = f(z);
w = zeros(1,41);
for i=1:length(z)
 w(i) = newton_interp(x,y,z(i)); % Interpolating using equidistance nodes
end
subplot(3,1,1)
plot(z,ynew,'ro',z,w,'b-.'); 
legend('y','w')

%% Interpolating using Chebyshev nodes
for i=1:21
 x(i) = 5 * cos((2*i-1)*pi/42); 
end
y = f(x);
tic
for i = 1:length(z)
 w(i) = newton_interp(x,y,z(i));
end
toc
y = f(z);
subplot(3,1,2)
plot(z,y,'ro',z,w,'b-')
legend('y','w-Chebyshev')

x = linspace(-5,5,21);
f = @(x) 1./(1+x.^2);
y = f(x);
tic
[ls,S] = polyfit(x,y,17);
toc 
d = ls(1) * x + ls(2) - y;
J = sum(d.^2);
v = polyval(ls,z);

subplot(3,1,3)
plot(z,w,'b-',z,v,'ro')
legend('w-Chebyshev','linear least square')







