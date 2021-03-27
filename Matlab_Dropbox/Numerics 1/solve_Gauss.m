function x=solve_Gauss(A,b)

%script to solve a linear system by using Gaussian elimination

% Check whether sizes of the matrices fit
[m,n]=size(A);
[p,q]=size(b);

if m ~= n
    error('The matrix A is not square')
end

if n~=p
    error('The number of rows in b does not match with the one in A')
end

[L,U]=LU_dec(A)

disp('Norm-error of the LU-decomposition (without pivoting) is')
norm(A-L*U)

b = inv(L) * b;

disp('Solution of the system is')
% x = inv(U) * b
% Using the solver for upper triangular system
x = solve_Ub(U,b);



