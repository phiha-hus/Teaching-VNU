% Ham so chung cho 3 phuong phap Hinh thang hien; Heun; Ralston
% [T,Y] = RK2(f,t0,tf,y0,h,c)
% c la vector he so 

% Hinh thang hien: c = [1/2 1/2 0 1]
% Heun
% Ralston 
function [T,Y] = RK2(f,t0,tf,y0,h,c)

c0 = c(1); c1 = c(2); p = c(3); q = c(4);

N = ceil((tf-t0)/h); T = linspace(t0,tf,N);
Y = [y0]; 

for i = 2:N
    ti = T(i); yi = Y(end);
    k0 = h .* f(ti,yi);
    k1 = h .* f(ti + p*h,yi + q*k0);
    yi_new = yi + c0 * k0 + c1 * k1;
    Y = [Y yi_new];
end
