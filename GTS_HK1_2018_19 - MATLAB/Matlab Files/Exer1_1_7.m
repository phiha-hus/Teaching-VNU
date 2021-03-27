% Exercise 7 (section1.1)
clear all; close all; clc

tol1 = [];
tol2 = [];
N = 20;
I = linspace(2,N,N);

for i=1:N
    %x = rand/2.^i;    
    x = 1/2.^i;    
    tol1 = [tol1 (1+x).^2-1];
    tol2 = [tol2 x.^2+2*x]; 
end

subplot(2,2,1)
semilogy(I,tol1,'r*-')    
title('Absolute error (1+x)^2-1')    

subplot(2,2,2)
semilogy(I,tol2,'b*-')    
title('Absolute error x^2+2*x')    

subplot(2,2,3)
plot(I,tol1-tol2,'r*-')
title('Difference tol1-tol2')    

subplot(2,2,4)
semilogy(I,tol2-tol1,'b*-')
title('Difference tol2-tol1')    