T = 2;
q0 = 120;
qT = 60;
q0p = 0;
qTp = 0;

A = [1, 0, 0, 0;
     1, T, T^2, T^3;
     0, 1, 0, 0;
     0, 1, 2*T, 3*T^2];
b = [q0; qT; q0p; qTp];

x = A \ b;

a0 = x(1);
a1 = x(2);
a2 = x(3);
a3 = x(4);

fprintf('a0 = %f, a1 = %f, a2 = %f, a3 = %f\n', a0, a1, a2, a3);

q = a0 + a1*t + a2*t.^2 + a3*t.^3;
f = t.^3 - 3*t.^2 + 5.*t.*sin((pi.*t - 5*pi)/4) + 3;

t = (0:200)*0.01;
plot(t, q);
hold on;
plot(t, f);
xlabel('T');
ylabel('q(T)');
title('Graph');
legend('q(x)', 'f(x)');

