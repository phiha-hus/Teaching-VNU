h = 1e-2;
t = 0:h:2;
y = zeros(size(t)); 
y(1) = 1;

f = @(t, y) (t.*y + exp(t).*(3.*t.^2 + 2.*t + 1))./(1+t);

tspan = [0, 2];

for i = 1:length(t)-1
    y_temp = y(i) + h*f(t(i), y(i));
    y(i+1) = y(i) + h/2*(f(t(i), y(i)) + f(t(i+1), y_temp));
end

y0 = 1;

%ode45
[t1,y1] = ode45(f, tspan, y0);

%ode15s
[t2,y2] = ode15s(f, tspan, y0);

subplot(1,3,1);
plot(t, y, '-o');
xlabel('t');
ylabel('y');
title('Heun');

subplot(1,3,2);
plot(t1,y1,'-o');
xlabel('t');
ylabel('y');
title('ode45');

subplot(1,3,3);
plot(t2,y2,'-o');
xlabel('t');
ylabel('y');
title('ode15s');
