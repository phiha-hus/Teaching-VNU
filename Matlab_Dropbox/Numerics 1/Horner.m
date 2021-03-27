function [y,time1,time2] = Horner(a,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% This function computes the value of the polynomial 
% P(x) = a0 +a1 * x + ... + an *x^n
% at the value x=r
% using the Horner scheme
% Written by Phi, Dec. 2016 
% Use the book Kincaid-Cheney or Thomas-Fink
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all; close all; clc

%a = input('Input the coefficients a0,a1,...,an of P correspond to the power 1,x,...,x^n as one row vector: ')
%r = input('Input the value c that you want to evaluate the polynomial P(c): ')

[~,n] = size(a)
b = zeros(1,n-1);

b(n-1) = a(n);

% First way, using one for loop
tic
for j=n-2:-1:1
    b(j) = b(j+1) * r + a(j+1);
end

format short
disp('The quotient is')
b
disp('Evaluation of p at x=r')
y = b(1) * r + a(1)
time1 = toc 

% Second way, use matrix equation
% [b1 ... b_{n-1}] = [b2 ... b_{n-1}] * r + [a2 ... an]

A = eye(n-1) - diag(r*ones(n-2,1),1) % Careful with diag(r*ones(n-2),1) ==> sizes do not match
size(A)
size(a(:,2:n))

 tic
 disp('The quotient is')
 b = A\a(:,2:n)'     % Use matrix division \ instead of inv(A) in order to avoid computing inv(A)
 disp('Evaluation of p at x=r')
 y = b(1) * r + a(1)
 time2 = toc

end


%%







