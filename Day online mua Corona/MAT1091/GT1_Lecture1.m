% Slide Giai Tich 1
% Minh hoa gioi han cua day so va cua ham so

x = 0:2:99;

y = 1./x;


figure(1); clf;
subplot(1,2,1)
semilogy(x,y,'rs')
grid on
xlabel('Chi so n cua u_n')
ylabel('Gia tri cua u_n')
legend('u_n')
title('Gioi han cua day so u_n = 1/n khi n -> vo cung. Viet lim u_n = 0')


y = x.^2+1;
subplot(1,2,2)
semilogy(x,y,'rs')
grid on
xlabel('Chi so n cua v_n')
ylabel('Gia tri cua v_n')
legend('v_n')
title('Gioi han cua day so v_n = n^2+1 khi n -> vo cung')




x = 0:2:100;
figure(2); clf;

subplot(1,2,1)
y = sin(x);
plot(x,y,'rs')
grid on
legend('sin(n)')
title('Do thi ham so sin(n)')

subplot(1,2,2)
y = cos(x);
plot(x,y,'rs')
grid on
legend('sin(n)')
title('Do thi ham so cos(n)')

