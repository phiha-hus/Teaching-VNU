pkg load control
pkg load signal

A = [-1 0 -1; 0 0 1; 1 -1 0]; B = [1;0;0]; C = [1 -1 0]; D = 0;
[N1,D1] = ss2tf(A,B,C,D,1);
[Z1,P1,K1] = tf2zp(N1,D1)
x = [1;1;1];

g = ss(A,B,C,D);
initial(g,x);