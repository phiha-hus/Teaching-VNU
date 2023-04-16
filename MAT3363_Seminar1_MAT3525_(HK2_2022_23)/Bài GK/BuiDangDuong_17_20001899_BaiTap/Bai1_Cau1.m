% Điều kiện ban đầu
T = 2; 
q0 = 120; 
qT = 60; 
q0_dot = 0;
qT_dot = 0; 

% Giải hệ phương trình để tìm các hệ số a0, a1, a2, a3
A = [1, 0, 0, 0; 1, T, T^2, T^3; 0, 1, 0, 0; 0, 1, 2*T, 3*T^2];
B = [q0; qT; q0_dot; qT_dot];
coefficients = linsolve(A, B); % Giải hệ phương trình tuyến tính

% Lấy các hệ số a0, a1, a2, a3 từ kết quả
a0 = coefficients(1);
a1 = coefficients(2);
a2 = coefficients(3);
a3 = coefficients(4);

% Hiển thị kết quả
disp("Hệ số của đa thức bậc 3 là:");
disp("a0 = " + a0);
disp("a1 = " + a1);
disp("a2 = " + a2);
disp("a3 = " + a3);
