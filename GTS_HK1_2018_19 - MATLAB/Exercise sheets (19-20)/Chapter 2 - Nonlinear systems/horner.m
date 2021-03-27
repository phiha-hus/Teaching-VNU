function [y,b] = horner(a,z)
%HORNER Horner algorithm
%   Y=HORNER(A,Z) computes
%   Y = A(1)*Z^N + A(2)*Z^(N-1) + ... + A(N)*Z + A(N+1)
%   using Horner's synthetic division algorithm.
n = length(a)-1;
b = zeros(n+1,1);
b(1) = a(1);
for j=2:n+1
   b(j) = a(j)+b(j-1)*z;
end
y = b(n+1);
b = b(1:end-1);
return
