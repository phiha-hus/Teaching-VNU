function test_fix_iter(index)

%format long
switch index

case 1 % Method converges
  f = @(x) exp(-x)
  x0 = 0.5  
case 2 % Method fails to converge
  f = @(x) 1 + 1/tan(x)
  x0 = 0.5  
case 3 % Converges fast - this is Newton method - quadratic convergence
  f = @(x) 1/2 * (x + 5/x)
  x0 = 2.5  
case 4
  f = @(x) (15*x.^2 - 24 * x + 13)/(4*x)
  alpha = 1; x0 = alpha + 1e-1
case 5  
  f = @(x) 3/4 * x + 1/(x.^3)
  alpha = sqrt(2);  
  x0 = alpha + 1e-1  
case 6  
  f = @(x) 27/sqrt(x) ; 
  x0 = 2.5  ; 
case 7
  f = @(x) 32 * x.^6 - 48 * x.^4 - + 18 * x.^2 - 1;
  x0 = 1;
end

% format long
%[x,n] = fix_iter(f,x0,100,1e-4)

%tic
%[x,n] = fix_iter(f,x0,100,1e-4)
%toc

%tic
[x,n] = fix_iter(f,x0,100,1e-4)
%toc