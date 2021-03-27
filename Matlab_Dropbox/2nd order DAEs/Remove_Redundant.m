% This function is written to remove hidden redundancy in matrix pair (P,Q)
% based on Lemma 2.12 (PhD Thesis Phi Ha)
% Old schemes in HaM12 and HaMS13 is not yet fully implemented.
% rq: rankQ
% S and Z1 will be determined based on the SVD of P2
% Copyright: @Phi Ha

function [S,Z1,Z2] = Remove_Redundant(P,Q)

[m1,n1]=size(P);
[m2,n2]=size(Q);
if n1~=n2
    error('Column sizes of P and Q do not match')
end

if (rank(Q) ~= 0)    
    [Uq,Sq,Vq] = svd(Q);

    rq = rank(Sq) ;
    %tol = max(size(Q))*eps(max(Sq));
    %tol = 1e-4;
    %rq = nnz(diag(Sq))
    U11 = Uq(1:end,1:rq);

    V11 = Vq(1:rq,1:end);
    V12 = Vq(rq+1:end,1:end);

    Pnew = P * Vq;
    P1 = Pnew(1:end,1:rq);
    P2 = Pnew(1:end,rq+1:end);

    [Up,Sp2,~]=svd(P2);
    rp2 = rank(Sp2);
    %tol = max(size(Q))*eps(max(Sq));
    %tol = 1e-4;
    %rp2 = sum(Sp2 > tol)
    
    S = Up(1:rp2,1:end);
    Z1 = Up(rp2+1:end,1:end);
    P11 = S * P1;
    P21 = Z1 * P1;
    P12 = S * P2;
    Z2 = - P21 * (Sq(1:rq,1:rq) \ U11');  % can we improve inv here?
else
   S = eye(m1); Z1 = zeros(0,m1); Z2 = zeros(0,m2);     
end

