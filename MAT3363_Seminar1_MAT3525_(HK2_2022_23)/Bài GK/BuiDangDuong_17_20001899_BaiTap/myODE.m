function dydt = myODE(t, y)
    dydt = (t.*y + exp(t).*(3.*t.^2 + 2.*t + 1))./(1 + t);
end
h = 1e-2; % Bước thời gian
tspan = 0:h:2; % Khoảng thời gian từ 0 đến 2 với bước h
y0 = 1; % Điều kiện ban đầu y(0) = 1

% Sử dụng phương pháp Heun
[t1, y1] = heun(@myODE, tspan, y0, h);
% Sử dụng ode45
[t2, y2] = ode45(@myODE, tspan, y0);
% Sử dụng ode15s
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t3, y3] = ode15s(@myODE, tspan, y0, options);
% Tính nghiệm chính xác y(t) = exp(t)*(t^2 + 1)
y_exact = exp(tspan).*(tspan.^2 + 1);

% Tính sai số tuyệt đối
error1 = abs(y1 - y_exact');
error2 = abs(y2 - y_exact');
error3 = abs(y3 - y_exact');

% Vẽ đồ thị sai số tuyệt đối
plot(tspan, error1, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Heun');
hold on;
plot(tspan, error2, 'r-', 'LineWidth', 1.5, 'DisplayName', 'ode45');
plot(tspan, error3, 'g-', 'LineWidth', 1.5, 'DisplayName', 'ode15s');
xlabel('t');
ylabel('Absolute Error');
title('Absolute Error vs. t');
legend;
grid on;
