f = @(x, y) x^2 * (x + y);

% Định nghĩa giới hạn tích phân
x1 = 0;
x2 = 2;
y1 = @(x) x;
y2 = 3;

% Tính tích phân kép
result = integral2(f, x1, x2, y1, y2);
