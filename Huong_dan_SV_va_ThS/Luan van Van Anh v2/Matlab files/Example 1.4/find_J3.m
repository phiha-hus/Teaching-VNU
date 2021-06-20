function [K3] = find_K(h, A, Ad)
global	h	A	Ad
tic

h = 0.1;	kp = 0.4451; K = 1.53; ki = 2.3046; tau = 0.0254;
A = [ 0 1; 0 -1/tau];
Ad = [0 0; -K * ki/tau -K * kp/tau]
a = A;
N = 1;
I = eye(2)
for	k = -N : N
k 
W0 = lambertw_matrix(0,-h*Ad*expm(h*A));
W = fsolve(@witer,W0);
S0 = W0/h	-	A;
lambda = eig(S0);
%pause
end
alpha = max(lambda); 

h1 = 1e-2;
t = 0 : h1 : h;
syms s
 J3 = int(norm(expm((-A-alpha*I)*t +A*h)*Ad), s, 0, t)

 K3 = max(J3)
toc



