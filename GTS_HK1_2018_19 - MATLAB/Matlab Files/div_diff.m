function v = div_diff(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function div_diff is an implementation of divided difference scheme
%                   in order to compute Newton interpolation polynomial
% Form: v = div_diff(x,y)
% x: interpolation nodes 
% y: value of the function at nodes in vector x
% v: Output vector, contains the coefficients of the Newton interpolation polynomial
% Copyright: Phi Ha, October 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(x) ~= length(y)
  error('Length of nodes mismatch')
end

n = length(x);  
A = zeros(n,n);
A(:,1) = y';

for j = 2:n
  for i = j:n
    diff = x(i) - x(i-j+1);
    A(i,j) = (A(i,j-1)-A(i-1,j-1))/diff ;
  end
end 

v = diag(A)';