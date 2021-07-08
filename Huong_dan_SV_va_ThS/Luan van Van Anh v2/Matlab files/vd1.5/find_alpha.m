function [alpha] = find_alpha(A,Ad,h)
%FIND_ALPHA Summary of this function goes here
%   Detailed explanation goes here

[n,~] = size(A);
m = n - rank(Ad);

Eig_Sk = [];  % Long vector Eig_Sk store all eigenvalues of Sk

for k = -m:1:m

k    
temp = -h*Ad*expm(h*A) ;
D_init = lambertw_matrix(k,temp)
% W = fsolve(@witer,W0,optimoptions('fsolve','Display','iter'));

options = optimoptions('fsolve','Display','off');
[Dk,fval,exitflag,output] = fsolve(@witer,D_init,options) ;


Sk = Dk/h - A  
Eig_Sk = [Eig_Sk eig(Sk)]; 

end

alpha = max(real(Eig_Sk)) ; 

return

