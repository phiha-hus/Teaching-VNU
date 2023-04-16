function [T,Y] = euler_an(f,t0,tf,y0,h)

if nargin == 4
    h = 1e-2;
end

N = ceil((tf-t0)/h); T = linspace(t0,tf,N);
Y = [y0]; 

for i = 2:N
    ti = T(i);
    yi = Y(end);
    g = @(y) y - yi - h * f(ti,y); % y_{i+1} = y_i + h * f(t_{i+1},y_{i+1})
    yi_new = fsolve(g,yi);
    Y = [Y yi_new];
end


