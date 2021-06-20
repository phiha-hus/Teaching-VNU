% Function find_K4 in 1D
% Equation dx(t)/dt = A x(t) + Ad x(t-h) 

function [K4] = find_K4(h,A,Ad,alpha)

a = -A; ad = -Ad;
N = 50 ; tf = 2 * h ; 

tol = 1e-6; 
K4 = J4(N,tf,h,a,ad,alpha); 
K4_new = J4(N,tf,h,a,ad,alpha);

while abs(K4-K4_new)>tol
    K4 = K4_new;
    N = 2 * N;
    K4_new = J4(N,tf,h,a,ad,alpha);
end
return


function J4N = J4(N,tf,h,a,ad,alpha)   
% This function find sup of J4(N,t) w.r.t t on [h,tf]
% N is fixed
   h1 = 1e-2; T = h : h1 : tf;    
   value = 0;       
   for k = -N:1:N
    Sk = lambertw_matrix(k,-ad*h*exp(a*h)) - a ;
    num = ad * exp(-Sk * h) * exp((Sk - alpha) * T) ; 
    den  = (1 - ad * h * expm(-Sk * h)) ; 
    value = value + num/den ; 
   end 
   J4N = max(abs(value)) ;   
return
