function [L,U]=LU_dec(A)
% Finding the LU-decomposition of A
% without pivoting
promt = 'Choose 1 for LU_decomposition with pivoting, 0 for without pivoting  ';
r = input(promt);

% Check whether sizes of the matrices fit
[n,m]=size(A);
if m ~= n
    disp('Becareful! The matrix is not square')
end

% Set tolerance to verify whether an element on the main diagonal of A is
% zero
tol=1e-6;

L = eye(n);
R = A;

if r==0
    disp('LU decomposition withput pivoting')

for k=1:n-1        
        if norm(R(k,k)) < tol
           error(strcat('The ',num2str(k),' elements on the main diagonal are nealy zero, i.e., smaller than the tolerance')) 
        else
           disp(strcat('For k= ',num2str(k)) ) ;          
           L((k+1):n,k) = R((k+1):n,k)/R(k,k);
           
           % Eliminate x_k in all the equations from (k+1)_th to n_th
           % Notice that L((k+1):n,k) is a column, and R(k,(k+1):n) is a row
           % so the product means that scale the row with factors in the
           % column (Tensor product)
           R((k+1):n,(k+1):m) = R((k+1):n,(k+1):m) - L((k+1):n,k) * R(k,(k+1):m);
        end
%         disp('The matrix is')
%         R
end
  
else
   disp('LU decomposition with pivoting')

for k=1:n-1        
        if norm(R(k,k:n)) < tol
           error(strcat('All the elements on the', num2str(k), 'columns are nealy zero, i.e., smaller than the tolerance')) 
        else
           % Finding the maximum element i in the k_th column for permuting 2 rows
           % and then permute the k_th row with the i_th row
           [~,I] = max(R(k,k:n));
           S = R(I,:);
           R(I,:) = R(k,:);
           R(k,:) = S;          
           
           disp(strcat('For k= ',num2str(k)) ) ;          
           L((k+1):n,k) = R((k+1):n,k)/R(k,k);
           
           % Eliminate x_k in all the equations from (k+1)_th to n_th
           % Notice that L((k+1):n,k) is a column, and R(k,(k+1):n) is a row
           % so the product means that scale the row with factors in the
           % column (Tensor product)
           R((k+1):n,(k+1):m) = R((k+1):n,(k+1):m) - L((k+1):n,k) * R(k,(k+1):m);
        end
%         disp('The matrix is')
%         R
end 

end
% Notice that only the uppertriangular part of A is used. So R is not in
% the upper triangular form. Thus, we need to extract it from R
U = triu(R);