f = @(x, y) x^2 * sin(y);

% Giới hạn tích phân cho biến x và y
x1 = 0;
x2 = 4;
y1 = 0;
y2 = pi;

% Tính tích phân kép
result = integral2(f, x1, x2, y1, y2);
