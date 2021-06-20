function [K2] = find_S_VA(h, A, Ad)
global	h	A	Ad

tic 
h = 1;	A = 1;	Ad = 1;	N = 20; t0 = 0 ; 
a = A; ad = Ad; 

W0 = lambertw_matrix(0,-h*Ad*expm(h*A));
alpha = real(W0/h - A)

S = [];
% T chay tren [h,\infty)
h1 = 1e-2;
T = h : h1 : 20*h; 
[~,nT] = size(T);

J2N = zeros(2*N+1);

for	k = -N : N  
    k; 
    
    W = fsolve(@witer,W0);
    Sk = (1/h) * W	+	A;
    lambda = eig(Sk);
    S = [S Sk];
    %pause    
    
    tuso_k = exp((Sk - alpha) * T)  ;
    mauso_k = 1 - ad * h * expm(-Sk * h) ;  
    phanso_k = tuso_k./mauso_k;
  
    tong_phan_so = 0;
    tong_phan_so = tong_phan_so + phanso_k  
   
   J2N(i) = J2Nt ;
end

K2 = max(J2N) ; 

toc
