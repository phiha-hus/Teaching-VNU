% Test the bisection function written by Phi
function test_newton(index)
  
switch index

case 1        
    f = @(x) exp(x)- x - 10;
    df = @(x) exp(x) - 1; 
    x0 = 1;
    %[c,err] = bisection(f,0,5,100,1e-6)
    %[c,err,~] = bisection(f,0,5,100,1e-6)    
case 2
    f = @(x) x - exp(-x);
    df = @(x) 1 + exp(-x);
    x0  = 1 ;
case 3
    f = @(x) tan(x)-x;
    df = @(x) 1/cos(x).^2 - 1;
    x0 = 1 ;
case 4
    f = @(x) sin(x)+1;
    df = @(x) cos(x);
    x0 = 1 ;
case 5
    f = @(x) sin(x) - x;
    df = @(x) cos(x) - 1;
    x0 = 1e-4 ;
case 6
    f = @(x) x.^6 - x - 1;
    df = @(x) 6 * x.^5 - 1;
    x0 = 1.5;    
case 7
    f  = @(x) [x(1)+x(2)+x(3)-3; x(1).^2+x(2).^2+x(3).^2-5; exp(x(1))+ x(1)* x(2) - x(1) * x(3)-1];
    df = @(x) [1 1 1; 2*x(1) 2*x(2) 2*x(3); exp(x(1))+ x(2) - x(3) x(1) -x(1)];
    % x0 = [0.1 1.2 2.5]';        
    x0 = [1 0 1]';
case 8
    f = @(x) x.^3 - (3+1e-16) * x.^2 + 3 * x -(1+1e-16);
    df = @(x) 3 * x.^2 - 2 * (3+1e-12) * x + 3;  
     x0 = 1; 
end

tic
x  = bisection(f,0,5,100,1e-6)
%[x,~,n]  = bisection(f,0,5)
toc 

tic  
[x,n]  = newton(f,df,x0)
toc 
  
%% Compare with the built in function fzero of Matlab.
%tic 
%% x_approx = fzero(f,x0)   % Can only solve scalar equation but not system
%x_approx = fsolve(f,x0) 
%toc

%disp("Absolute error in comparison with fzero function in Matlab/Octave is: ")         
%abs(x-x_approx)