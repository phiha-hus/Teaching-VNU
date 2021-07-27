function [K1, alpha] = find_K(h, K, kp, ki, tau)
global	h	A	Ad

h = 0.1; kp = 0.4451; K = 1.53; 
ki = 2.3046; tau = 0.0254;

A = [ 0 1; 0 -1/tau];
Ad = [0 0; -K * ki/tau -K * kp/tau]

n = 2;
N = 1; 
for	k = -N : N
k 
W0 = lambertw_matrix(0,-h*Ad*expm(h*A));
W = fsolve(@witer,W0);
S0 = W0/h	-	A;
lambda = eig(S0);
%pause
end
alpha = max(lambda) 

h1 = 1e-2;
t = 0 : h1 : h;
J1 = norm(expm((-A-alpha*eye(n))*t))
K1 = max (J1)
   
  
end


    
    


