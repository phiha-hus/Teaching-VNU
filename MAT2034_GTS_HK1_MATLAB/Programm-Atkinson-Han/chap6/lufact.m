function [x, lu]=lufact(A,b)

% This will do an LU factorization of A and then call solve_tri
% to solve the linear system Ax=b.  Row pivoting is used, so that
% the eventual LU factorization returned in lu has been affected
% by row permutations.  This is meant to illustrate Gaussian 
% elimination and to illustrate the use of Matlab, possibly in
% not the most efficient manner.  The student should feel free to
% improve on this code.

n=size(A,1);

for k=1:n-1
   % Find the maximal element in the pivot column, below the 
   % pivot position, along with the index of that maximal element.
   
   [col_max   index] = max(abs(A(k:n,k)));
   index = index + k-1;
   
   if index ~= k
      % Switch rows k and index, in columns k thru n.  Do similarly
      % for the right-hand side b.
      
      tempA = A(k,k:n);
      A(k,k:n) = A(index,k:n);
      A(index,k:n) = tempA;
      
      tempb = b(k);
      b(k) = b(index);
      b(index) = tempb;
   end
   
   % Form the needed multipliers and store them into the pivot
   % column, below the diagonal.
   A(k+1:n,k) = A(k+1:n,k)/A(k,k);
   
   % Carry out the elimination step, first modifying the matrix,
   % and then modifying the right-hand side.
   for i=k+1:n
      A(i,k+1:n) = A(i,k+1:n) - A(i,k)*A(k,k+1:n);
   end
   b(k+1:n) = b(k+1:n) - A(k+1:n,k)*b(k);
end

% Solve the upper triangular linear system.
x = solve_tri(A,b);
lu = A;
