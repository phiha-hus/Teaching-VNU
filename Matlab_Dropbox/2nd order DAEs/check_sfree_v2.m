function [is_sfree,MNP] = check_sfree_v2(A2,A1,A0)

[~,n] = size(A2);
[r2,r1,~,A2,A1,A0] = BlkUpTri(A2,A1,A0);

% Select the set of 1st order scalar difference equations in the strangeness-free form
% Similar to the process of Removing hidden redundancy in A12 and A03
[S,~,~] = Remove_Redundant(A1(r2+1:end,:),A0(r2+r1+1:end,:));

% Select the set of 1st order scalar difference equations in the strangeness-free form
% Similar to the process of Removing hidden redundancy in [A12;A03] and A21
A_12_03 = [S*A1(r2+1:end,:) ; A0(r2+r1+1:end,:)];
[U,~,~] = Remove_Redundant(A2(1:r2,:),A_12_03);

MNP = [U*A2(1:r2,:); A_12_03];
is_sfree = (rank(MNP) == n);    

if is_sfree==0
    disp('bigger strangeness-index')
else
    disp('Strangeness-index STOPS here')
end

