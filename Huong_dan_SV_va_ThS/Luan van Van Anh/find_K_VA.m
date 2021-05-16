function [K4] = find_K_VA(h, A, Ad)
global	h	A	Ad
tic

h=1;	A = 1;	Ad = 1;	N = 10; t0 = 0
a = A; ad = Ad; 

W0 = lambertw_matrix(0,-h*Ad*expm(h*A));
alpha = real(W0/h - A);

S = [];
% T chay tren [h,\infty)
h1 = 1e-1;
T = h : h1 : 20*h; 

J4N = zeros(size(T));

for i = 1:length(T)
    t = T(i);
  J4Nt = 0; 
  for	k = -N : N  
    k; 
    W = fsolve(@witer,W0);
    Sk = (1/h) * W	+	A;
    lambda = eig(Sk);
    S = [S Sk];
    %pause
    
 
    mauso_k = 1 - ad * h * expm(-Sk * h) ;
    tuso_k = ad * exp((Sk - alpha) * t)  ;
    phanso_k = tuso_k/mauso_k;
  
    syms s
    I_k = int(phanso_k * expm (-Sk * s), s , 0 , h)

    tong_Ik = 0;
    tong_Ik = tong_Ik + I_k  ;
   end
   J4Nt = norm(tong_Ik) ; 
   
   J4N(i) = J4Nt ;
end

K4 = max(J4N) ; 
toc



