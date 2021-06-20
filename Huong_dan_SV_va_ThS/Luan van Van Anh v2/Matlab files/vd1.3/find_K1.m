% Function find_K1 in 1D
function [K1] = find_K1(h,A,alpha)

a = A;

K1 = 1 ;  % if a <= alpha
if a > alpha 
    K1 = expm(a-alpha)*h ;
end

return