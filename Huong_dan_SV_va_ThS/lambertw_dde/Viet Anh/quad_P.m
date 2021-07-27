function [N,value] = quad_P(f,t0,tf,tol)

if nargin == 3
     tol = 1e-6;
end

N = 1e+3; 

tic

value_1 = quad_P0(f,t0,tf,N) ;
value_2 = quad_P0(f,t0,tf,2 * N);

while abs(value_1 - value_2) > tol
    N = 2 * N ;     
    value_1 = value_2 ;
    % 2nd way - I thought it is faster - but not - why?
    % h = (tf-t0)/N ; 
    % temp = quad_P0(f,t0+h,tf-h,N/2) ;  
    % value_2 = value_1/2 + h .* temp ;
    
    % 1st way, simpler
    value_2 = quad_P0(f,t0,tf,2*N) ;  
end

value = value_2 ; 

toc
return



function value = quad_P0(f,t0,tf,N)
    h = (tf-t0)/N;

    T = t0 + h * [0:N];

    f_value = f(t0) ;

    for i = 1:N
        f_value = f_value + f(T(i));
    end

    value = f_value .* h ;