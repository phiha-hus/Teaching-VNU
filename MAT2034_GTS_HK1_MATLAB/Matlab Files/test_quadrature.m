function test_quadrature(index)

% The thing that I can observe here is that computing I_{2n} instead of I_{n+10}
% does not make the problem much faster.
% May be for TRAPEZOIDAL RULE then computing I_{2n} instead of I_{n+1}
% BUT for SIMPSON THAN IT IS NOT FASTER TO COMPUTE USING I_{2n} instead of 
% $I(n=100)$ ------> so may be INCREASE N BY 20 EVERYTIME 

% I GUESS THE BETTER SOLUTION IS ADAPTIVE_SIMPSON

switch index
    case 1
        f = @(x) x.^2; a = 1; b = 10;
    case 2
        f = @(x) 4./(1+x.^2); a = 0; b = 1;
    case 3
end

% Test the trapezoidal rule - self implement - vs. Matlab built in function
% trapz

output_precision(7)

disp('Test the performance of trapezoidal: ')
[value,~,n] = trapezoidal(f,a,b)

disp('Test the performance of built in function trapez: ')
X = linspace(a,b,2.^n-1); Y = f(X);
tic
trapz(X,Y)
toc

disp('Test the performance of simpson')
tic 
[value,err,n]=simpson(f,a,b)
toc 

disp('Test the performance of built in function simpson')
tic
quadv(f,a,b)
toc

disp('Test the performance of simple simpson')
tic
simpson_simple(f,a,b,100)
toc
