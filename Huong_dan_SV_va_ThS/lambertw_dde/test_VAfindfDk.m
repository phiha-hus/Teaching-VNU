clear all; close all; clc

%% khi test thi doi khi phai dung nhung vi du that de truoc
% Test xong thi thay code pp Newton minh viet la dung roi

D_init = [1.1 0; 0 1.9];
f = @(D_k) D_k.^2 - [1 0;0 4];
imax = 10000; tol  = 1e-6; 
[Dk,i] = VAfind_Dk(f, tol, imax, D_init) ;


%% VA viet pp Newton khong sai dau, co dieu vi Jacobian qua lon 
% nen phai dung phuong phap khac Newton thoi - 
% Dung fsolve trong Matlab la duoc

k = 0;
tau     =  1;                                
A     =  [-1 -3;2 -5]; 
Ad    =  [1.66 -0.697;0.93 -0.330];
imax = 1e+4; tol  = 1e-6; 

% Day la VA viet nay
D_init = lambertw(k, tau*Ad*expm(tau*A)) ;  
f = @(Dk) Dk * expm(Dk - A * tau) - Ad * tau ;
[Dk,i] = VAfind_Dk(f, tol, imax, D_init) ;

% Day la thay dung fsolve nhe
options = optimoptions('fminunc'); 
options.TolX = 1e-6;
options.MaxIter = 1000;
Dk = fsolve(f,D_init,options)


