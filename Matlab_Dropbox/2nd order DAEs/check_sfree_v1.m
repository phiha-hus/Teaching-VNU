% Verify whether a given SiDE contains its strangeness-free formulation or
% not
% If yes, then we select the strangeness-free formulation by using
% supporting functions BlkUpTri & Remove_Redundant
% Copyright: @Phi Ha

function [is_sfree,M,N,P] = check_sfree_v1(A2,A1,A0)

[~,n] = size(A2);
%disp('Run the block triangular formulation')
[r2,r1,r0,A2,A1,A0] = BlkUpTri(A2,A1,A0);

A21 = A2(1:r2,:);
A11 = A1(1:r2,:); A12 = A1(r2+1:r2+r1,:);
A01 = A0(1:r2,:); A02 = A0(r2+1:r2+r1,:); A03 = A0(r2+r1+1:r2+r1+r0,:);

is_sfree = (rank([A21; A12; A03]) == n);

if is_sfree==0
    M = []; N = []; P = [];
    disp('bigger strangeness-index')
elseif is_sfree==1    
    disp('Now we need to select the strangeness-free form')
    % The fact that really suprise me is that we cannot use the
    % Remove_Redundant above to select the s-free-form
    
    % Select the set of 1st order scalar difference equations in the strangeness-free form
    % Removing hidden redundancy in A12 and A03
    [U,~,~] = Remove_Redundant(A12,A03);
    [d1,~] = size(U);
    A12_new = U * A12;
    
    % Removing hidden redundancy in A21 and [A12;A03]
    [S,~,~] = Remove_Redundant(A21,[U * A12; A03]);
    d2 = size(S,1);
    A21_new = S*A21;

    % Construct the strangeness-free formulation
    M = [S * A21; zeros(n-d2,n)];
    N = [S * A11; U * A12; zeros(n-d1-d2,n)];
    P = [S * A01; U * A02; A03];
end

