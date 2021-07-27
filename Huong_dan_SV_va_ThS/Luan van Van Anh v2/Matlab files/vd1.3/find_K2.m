% Function find_K2 in 1D
% Equation dx(t)/dt = A x(t) + Ad x(t-h) 

function [K2] = find_K2(h,A,Ad,alpha)

a = -A; ad = -Ad;
N = 50 ; tf = 2 * h ; 

tol = 1e-6; 
K2 = J2(N,tf,h,a,ad,alpha); 
K2_new = J2(N,tf,h,a,ad,alpha);

while abs(K2-K2_new)>tol
    K2 = K2_new;
    N = 2 * N;
    K2_new = J2(N,tf,h,a,ad,alpha);
end
return


function J2N = J2(N,tf,h,a,ad,alpha)   
% This function find sup of J2(N,t) w.r.t t on [h,tf]
% N is fixed
   h1 = 1e-2; T = h : h1 : tf;    
   value = 0;       
   for k = -N:1:N
    Sk = lambertw_matrix(k,-ad*h*exp(a*h)) - a ;
    num = exp((Sk - alpha) * T) ; 
    den  = (1 - ad * h * expm(-Sk * h)) ; 
    value = value + num/den ; 
   end 
   J2N = max(abs(value)) ;   
return


