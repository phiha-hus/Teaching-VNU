
%a Ham truyen
A = [-1/6 0 -1/3; 0 0 1; 0.5 -0.5 -0.5];
B = [1/6 1/3; 0 0; 0 0];
C = [1 -1 -1; -0.5 0 0];
D = [0 0; 0.5 0];
sys = ss(A, B, C, D);
sys_as_tf = tf(sys)
%a Khong diem cua he
z = tzero(sys_as_tf)
%a Cac cuc
P = pole(sys_as_tf)
%b
G = ss(A, B, C, D)
x0 = [1; 1; 1]
initial(G, x0)
grid
%c
%step(sys)
%impulse(sys)
