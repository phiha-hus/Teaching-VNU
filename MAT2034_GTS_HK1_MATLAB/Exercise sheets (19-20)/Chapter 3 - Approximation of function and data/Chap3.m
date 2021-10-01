% Chapter 3: Quarteroni
% Interpolation and data approximation
% Part 1: 
% Approximation by Taylor expansion, 
% Approximation by Lagrange polynomial and Polynomial in least square sense
% Approximation with special nodes to avoid Runge phenomenon

clear all; close all; clc

%% Approximation by Taylor expansion
f = @(x)exp(x)
taylortool(f)
grid on

%% Approximation by Lagrange polynomial and Polynomial in least square sense

K = -0.67;
x = [65 35 5 -25 -55];
y = [-3.1 -3.32 -3.02 -3.2 -3.25];

n = 4;

format short e

% if n = length(x)-1 then it is Lagrange interpolation
% otherwise it is an approximated n-order polynomial in least square sense
c = polyfit(x,y,n)

% z = linspace(x(1),x(end),100)
% p = polyval(c,z)
% plot(z,p,'k-.',x,y,'rs','linewidth',3);
% grid on

x1 = 65:-10:-55;
y1 = [-3.1 -3.22 -3.22 -3.32 -3.17 -3.07 -3.02 -3.02 -3.12 -3.2 -3.35 -3.37 -3.25];
z = x1;
p = polyval(c,z);
plot(z,p,'k-.',x1,y1,'rs','linewidth',3);
grid on

%% Approximation with special nodes to avoid Runge phenomenon
clear all; close all; clc

a= -5; b= 5; 
f = @(x) 1./(1+x.^2);

err1 = []; err2 = []; err3 = [];
S = 10:10:100;

for n = 10:10:100;

    % Chebyshev-Gauss-Lobatto
    xc = (a+b)/2 - (b-a)/2 * cos(pi * [0:n]/n);
    yc = f(xc);
    c1 = polyfit(xc,yc,n);

    % Chebyshev-Gauss
    xg = (a+b)/2 - (b-a)/2 * cos(pi/2 * (2*[0:n]+1)/(n+1));
    yg = f(xg);
    c2 = polyfit(xg,yg,n);

    % Equidistance
    xe = linspace(a,b,n);
    ye = f(xe);
    c3 = polyfit(xe,ye,n);

    % Uoc luong sai so
    x = linspace(-5,5,1000);
    y = f(x);

    pc = polyval(c1,x);
    pg = polyval(c2,x);
    pe = polyval(c3,x);

    err1 = [err1 max(pc-y)];
    err2 = [err2 max(pg-y)];
    err3 = [err3 max(pe-y)];

end

semilogy(S,err1,S,err2,S,err3)
grid on
legend('C-G-L','C-G','equi.distance')






