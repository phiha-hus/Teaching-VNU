% Chapter 3: Quarteroni
% Interpolation and data approximation
% Part 2: 
% Trigonometric interpolation
% Piecewise linear and spline interpolation
% Least square interpolation

clear all; close all; clc
%% Trigonometric interpolation
% Example 3.6

a = 0; b = 2*pi;
n = 9; f =@(x)  x .* (x-2*pi) .* exp(-x);

x = (b-a) .* [0:n]/(n+1);
y = f(x);
z = interpft(y,100)

t = linspace(a,b,100);

plot(t,f(t),'bs',t,z)
legend('f(x)','interpft'); hold on
grid on

plot(x,f(x),'r*')

%% Piecewise linear interpolation
clear all; close all; clc

n = 1e+4;
x = rand(n,1); x = sort(x);
y = rand(n,1);
z = rand(100 * n,1);

tic
interp1(x,y,z);
toc

tic
interp1q(x,y,z);
toc

%% Spline interpolation
clear all; close all; clc

x = [ -55:10:65];
y = [ -3.25 -3.37 -3.35 -3.2 -3.12 -3.02 -3.02 ...
-3.07 -3.17 -3.32 -3.3 -3.22 -3.1];
zi = [ -55:1:65];

tic
s = spline(x,y,zi );  % Matlab built in function spline
toc

tic
s3 = cubicspline(x,y,zi);
toc

type = 0; der = [1 2];
tic
s3 = cubicspline(x,y,zi,type,der);
toc

%% Interpolation in 2D, 3D

[x,y]= meshgrid (0:0.2:1 ,0:0.2:1);
z= sin (2* pi*x ).* cos (2* pi*y );
plot3(x,y,z)
grid on

xi = [0:0.05:1]; yi =[0:0.05:1];
[xf ,yf ]= meshgrid (xi ,yi );
pi3= interp2 (x,y,z,xf ,yf );








