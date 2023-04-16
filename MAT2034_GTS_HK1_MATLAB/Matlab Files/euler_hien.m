function [T,Y] = euler_hien(f,t0,tf,y0,h)

if nargin == 4
    h = 1e-2;
end

N = ceil((tf-t0)/h); T = linspace(t0,tf,N);
Y = [y0]; 

for i = 2:N
    ti = T(i);
    yi = Y(end);
    yi_new = yi + h * f(ti,yi); % y_{i+1} = yi_new
    Y = [Y yi_new];
end


