% This function construct a strangeness-free formulation of
% Linear Time Invariant 2nd order Singular Difference Equations (SiDEs)
% Input: 3 matrices A, B,C are system's coefficients of SiDEs
% Output: ka (kappa), M, N,P where
% ka (kappa): the strangeness-index of the matrix triple
% M,N, P are matrices in the strangeness-free form, without the set of
% redundant equations
% Copyright: @Phi Ha

function [ka,M,N,P]=sfree_triple(A,B,C)

[m,n] = size(A); [m1,n1] = size(B); [m2,n2] = size(C);
% Check match sizes
if (n~=n1) or (m~=m1)
    error('Sizes do not match')
end
if (n~=n2) or (m~=m2)
    error('Sizes do not match')
end

disp('We start with the strangeness-index = 0')
ka = 0; 
[is_sfree,M,N,P] = check_sfree_v1(A,B,C);
%[is_sfree,M,N,P] = check_sfree_v2(A,B,C);

if is_sfree==1
    disp('The strangeness-index is 0')    
end
    
while (ka<4) && (is_sfree==0)
    disp('Increase the index by 1 as: ')   
    ka=ka+1
    
    % Form the inflated matrix
    if ka==1
        Y = [C B A zeros(m,n);zeros(m,n) C B A];
    elseif ka==2
        Y = [C B A zeros(m,2*n);zeros(m,n) C B A zeros(m,n);zeros(m,2*n) C B A];
    elseif ka==3
        Y = [C B A zeros(m,3*n); zeros(m,n) C B A zeros(m,2*n); zeros(m,2*n) C B A zeros(m,n); zeros(m,3*n) C B A];
    else 
        disp('EITHER THE STRANGENESS INDEX IS BIGGER THAN 3 OR THERE IS NOT UNIQUE SOLUTION.')
    end
    
    W = null(Y(:,3*n+1:end)');
    
    % Update P,N,M for checking whether is_sfree = 0
    P = W' * Y(:,1:n);                
    N = W' * Y(:,n+1:2*n);
    M = W' * Y(:,2*n+1:3*n);        
   [is_sfree,M,N,P] = check_sfree_v1(M,N,P);    
   %[is_sfree,M,N,P] = check_sfree_v2(M,N,P)    
end

if ka>0
disp('The strangeness-index is: ')
ka
end


%-------------------------------------------------------------------------%
% Supporting functions.
%-------------------------------------------------------------------------%
