% This function construct a strangeness-free formulation of
% Linear Time Varying 2nd order Singular Difference Equations (SiDEs)
% Input: 3 matrix-valued function A, B,C are system's coefficients of SiDEs
%        and the discret time point t=n
% Output: ka (kappa), M, N,P where
%         ka (kappa): the strangeness-index of the SiDE at the time t=n
%         M,N, P are matrices in the strangeness-free form at the time t=n, 
%         without the set of redundant equations
% Copyright: @Phi Ha

function [ka,M,N,P]=sfree_triple_LTV(n,A,B,C)

An = A(n); Bn = B(n); Cn=C(n);

[m0,d0] = size(An); [m1,d1] = size(Bn); [m2,d2] = size(Cn);
% Check match sizes
if (d0~=d1) or (m0~=m1)
    error('Sizes do not match')
end
if (d0~=d2) or (m0~=m2)
    error('Sizes do not match')
end

disp('We start with the strangeness-index = 0')
ka = 0; 
[is_sfree,M,N,P] = check_sfree_v1(An,Bn,Cn);
%[is_sfree,M,N,P] = check_sfree_v2(An,Bn,Cn);

if is_sfree==1
    disp('The strangeness-index is 0')    
end
    
while (ka<4) && (is_sfree==0)
    disp('Increase the index by 1 as: ')   
    ka=ka+1
    
    % Form the inflated matrix
    if ka==1
        Y = [Cn Bn An zeros(m0,d0);zeros(m0,d0) C(n+1) B(n+1) A(n+1)];
    elseif ka==2
        Y = [Cn Bn An zeros(m0,2*d0);zeros(m0,d0) C(n+1) B(n+1) A(n+1) zeros(m0,d0);zeros(m0,2*d0) C(n+2) B(n+2) A(n+2)];
    elseif ka==3
        Y = [Cn Bn An zeros(m0,3*d0);zeros(m0,d0) C(n+1) B(n+1) A(n+1) zeros(m0,2*d0);...
            zeros(m0,2*d0) C(n+2) B(n+2) A(n+2) zeros(m0,d0); zeros(m0,3*d0) C(n+3) B(n+3) A(n+3)];        
    else 
        disp('EITHER THE STRANGENESS INDEX IS BIGGER THAN 3 OR THERE IS NOT UNIQUE SOLUTION.')
    end
    
    W = null(Y(:,3*d0+1:end)');
    
    % Update P,N,M for checking whether is_sfree = 0
    P = W' * Y(:,1:d0);                
    N = W' * Y(:,d0+1:2*d0);
    M = W' * Y(:,2*d0+1:3*d0);        
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
