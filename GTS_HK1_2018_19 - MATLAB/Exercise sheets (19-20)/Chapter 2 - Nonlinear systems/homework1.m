% Bai lam cua SV Vu Manh Cuong
% MS 16002504

function homework1()
    clear all; close all; clc
    %format long
    
    %ex2(300, 3.5*10^-7, 0.401, 42.7*10^-6, 1000, 1.3806503*10^-23);
    %ex3(1,1,9.8)
    ex12(9.8, 1, 10);
    %ex13(6000, 1000, 5);
    %ex14(8,10,3*pi/5);
end

function [V] = ex2(T, p, a, b,N, k)
    f2 = @(V) (p + a*(N/V)^2)*(V-N*b) - k*N*T;    
    V = solve(f2,0)
end

function omega = ex3(s, t, g)
    f3 = @(w) g/(2*w^2)*(sinh(w*t) - sin(w*t)) - s;
    omega = solve(f3,0);
end

% not complete
function ex9(a1, a2, a3, a4)
    beta = 0;
    f9 = @(x) a1/a2*cos(beta) - a1/a4*cos(x) - cos(beta - x) 
    + (a1^2+a2^2-a3^2+a4^2)/(2*a2*a4)
end

function alpha = ex12(g,h,v0)
    f12 = @(alpha) sin(alpha) - sqrt(2*g*h/(v0^2));
    df12 = @(alpha) cos(alpha);
    alpha = newton(f12,df12,0, 10^-5, 1000)
end

function r = ex13(M, v , n, tolerance)
    if nargin < 4
        tolerance = 10^-12;
    end
    syms r;
    f13 = M - v.*(1.+r)./r.*((1.+r).^n -1)
    df13 = diff(f13)
   
    f13 = matlabFunction(f13)
    df13 = matlabFunction(df13)
    r = newton(f13, df13, -1, tolerance, 1000)
end

function alpha = ex14(l1, l2, gama)
    syms alpha
    f14 = l2*cos(pi-alpha-gama)/(sin(pi-gama-alpha))^2 - l1*cos(alpha)/(sin(alpha))^2;
    df14 = diff(f14);
    
    f14 = matlabFunction(f14);
    df14 = matlabFunction(df14);
    
    fplot(f14,-10:1:10)
    
    alpha = newton(f14,df14, -1, 10^-5, 1000);
end