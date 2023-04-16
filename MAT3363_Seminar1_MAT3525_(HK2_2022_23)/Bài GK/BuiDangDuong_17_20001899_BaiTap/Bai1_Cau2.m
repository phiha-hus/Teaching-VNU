
t = 0:0.01:2; 
q = 120 + 0*t - 45*t.^2 + 15*t; % Đa thức bậc 3 đã cho

n = 3; % Bậc của đa thức
p = polyfit(t, q, n); % Hàm polyfit để xác định các hệ số

% Vẽ đồ thị q(T)
plot(t, q, 'b', 'LineWidth', 2);
xlabel('T');
ylabel('q(T)');
title('Đồ thị q(T)');
grid on;

