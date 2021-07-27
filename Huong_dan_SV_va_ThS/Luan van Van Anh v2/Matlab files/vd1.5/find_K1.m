function [K1, alpha] = find_K1(h, K, kp, ki, tau)
global	h	A	Ad

W0 = lambertw_matrix(0,-h*Ad*expm(h*A));
W = fsolve(@witer,W0);
S0 = W0/h - A;
lambda = eig(S0);
alpha = max(lambda)

n = 2;              % Dimension of the system
N = n - rank(Ad);   % Dimension of nullity(Ad)
h1 = 1e-2; t = 0:h1:h; % Discretize time interval

% Norm of J1 if t = 0
K1 = 1;

for i = 2:length(t)
    J1 = norm(expm((-A-alpha*eye(n)) * t(i) ));
    K1 = max(K1,J1);
end
   
return

    
    


