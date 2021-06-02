% Ch?n R = L = C = 1
%a
A = [-1 1;-1 -1];
B = [0; 1];
C = [-1 1];
D = [0];
[N1, D1] = ss2tf(A, B, C, D, 1)
N = [0 1 0];
D = [1 2 2];
sys = tf(N, D)
[Z1, P1, K1] = tf2zp(N1, D1)

%b
A = [-1 1;-1 -1];
B = [0; 1];
C = [-1 1];
D = [0];
x0 = [1; 1];
G = ss(A, B, C, D);
initial(G, x0)
grid

%c
sys = tf(1, [1 2 2]);
step(sys)
impulse(sys)
