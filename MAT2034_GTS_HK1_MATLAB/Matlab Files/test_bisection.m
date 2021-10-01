% Test the bisection function written by Phi
function test_bisection(index)
  
switch index

case 1        
    f = @(x) exp(x)- x - 10;    
case 2
    f = @(x) x - exp(-x);
case 3
    f = @(x) tan(x)-x;
case 4
    f = @(x) sin(x)+1;
case 5
    f = @(x) (x.^2-2*x+1)/(x.^2+x-2);
case 6
    f = @(x) (x.^2-2*x+1)/(x.^2-x-2);    % Example in VHLinh's exercise sheet 
case 7
    f = @(x) x.^3 - 3*x.^2 + 3*x - 1;
    f2 = @(x) - 1 + 3*x - 3*x.^2 + x.^3;
    g = @(x) -1 + x*(3 + x * (-3+x));
end

% c = bisection(f,0,5,100,1e-6)
% [c,err] = bisection(f,0,5,100,1e-6)
% [c,err,~] = bisection(f,0,5,100,1e-6)

% [c,err,n] = bisection(f,0,3)

% Compare with the built in function fzero of Matlab.
% c2 = fzero(f,[0 3])

% Compare different way of writing the function f ==> Case 7
tic 
[c,err,n] = bisection(f,0,3,1e+2,1e-6);
toc 

tic 
[c,err,n] = bisection(f2,0,3,1e+2,1e-6);
toc 

tic 
[c,err,n] = bisection(g,0,3,1e+2,1e-6);
toc


 

