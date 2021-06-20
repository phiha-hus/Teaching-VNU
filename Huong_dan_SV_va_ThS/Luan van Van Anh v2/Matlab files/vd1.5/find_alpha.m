function [alpha] = find_alpha(h,A,Ad)
%FIND_ALPHA Summary of this function goes here
%   Detailed explanation goes here

W0 = real(lambertw_matrix(0,-h*Ad*expm(h*A)))
% W = fsolve(@witer,W0,optimoptions('fsolve','Display','iter'));

W = fsolve(@witer,W0) 

S0 = W./h - A;

alpha = max( eig(S0) );

return

