% Test example
%% Example 1: Euler hien vs. Euler an (stiff)
f = @(t,y) - 30 * y + 0 * t;
t0 = 0; tf = 5; h = 1e-1; y0 = 1;
x = @(t) exp(-100*t); % Nghiem chinh xac

[T1,Y1] = euler_hien(f,t0,tf,y0,h);
Err_1 = abs(Y1-x(T1));

[T2,Y2] = euler_an(f,t0,tf,y0,h);
Err_2 = abs(Y2-x(T2));

figure(1); clf;
subplot(1,2,1)
plot(T1,Err_1,'r-')
legend('Euler hien')
grid on

subplot(1,2,2)
plot(T2,Err_2,'b-')
legend('Euler an')
grid on

%% Example 2: Euler an vs. hinh thang hien (stiff)

f = @(t,y) - 30 * y + 0 * t;
t0 = 0; tf = 5; h = 1e-1; y0 = 1;
x = @(t) exp(-100*t); % Nghiem chinh xac

% c = [0 1 0.5 0.5] % Phuong phap Euler cai tien
c = [0.5 0.5 1 1]; % Phuong phap hinh thang hien
[T1,Y1] = RK2(f,t0,tf,y0,h,c);
Err_1 = abs(Y1-x(T1));

[T2,Y2] = euler_an(f,t0,tf,y0,h);
Err_2 = abs(Y2-x(T2));

figure(2); clf;
subplot(1,2,1)
plot(T1,Err_1,'r-')
legend('Hinh thang hien')
grid on

subplot(1,2,2)
plot(T2,Err_2,'b-')
legend('Euler an')
grid on

%% Example 3: Hinh thang hien vs. hinh thang an (stiff)

f = @(t,y) [(-1) * y(1) + 2 * y(2); -30 * y(2)];

t0 = 0; tf = 5; h = 1e-1; y0 = [1 exp(2)]';
x = @(t) expm([-1 2;0 -30]*t); % Nghiem chinh xac

% c = [0 1 0.5 0.5] % Phuong phap Euler cai tien
c = [0.5 0.5 1 1]; % Phuong phap hinh thang hien
[T1,Y1] = RK2(f,t0,tf,y0,h,c);
Err_1 = abs(Y1-x(T1));

[T2,Y2] = euler_an(f,t0,tf,y0,h);
Err_2 = abs(Y2-x(T2));

figure(3); clf;
subplot(1,2,1)
plot(T1,Err_1,'r-')
legend('Hinh thang hien')
grid on

subplot(1,2,2)
plot(T2,Err_2,'b-')
legend('Euler an')
grid on
