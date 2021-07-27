function [K1] = J1(h)
%J1 Summary of this function goes here
% Return K1, see eq (47) in Dua12 
%   Detailed explanation goes here

global A Ad alpha

h1 = 1e-6;
t = 0:h1:h;

[n,~] = size(A) ;

K1 = 1; 
for i = 2:length(t)
    temp = expm( (-A-alpha*eye(n)) * t(i) ) ;
    K1 = max( K1,norm(temp));
end

return

