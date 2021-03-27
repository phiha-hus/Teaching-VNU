function value = neville_interp(x,y,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neville Interpolation Using Divided Difference Scheme
% Form: value = neville_interp(x,y,w)
% x: interpolation nodes 
% y: value of the function at nodes in vector x
% w : point where we need to compute f(x) using Neville interpolation scheme 
% value: Output value of f(x)
% Copyright: Phi Ha, October 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(x) ~= length(y)
  error("Length of nodes mismatch")
end

n = length(x);
S = zeros(n,n);
S(:,1) = y';

for j = 2:n
  for i = j:n
    diff = x(i) - x(i-j+1);
    S(i,j) = ( (w-x(i-j+1)) * S(i,j-1) - (w-x(i)) * S(i-1,j-1) )/diff;           
  end  
end

value = S(n,n);

