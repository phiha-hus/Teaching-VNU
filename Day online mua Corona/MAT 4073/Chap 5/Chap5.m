% Chapter 5: Quarteroni - Linear systems - Giai gan dung he tuyen tinh Ax=b

%% Direct method : cac phuong phap truc tiep
% Problems 5.1 and 5.3 are very interesting
% 4 important decompositions: LU/PLU, Cholesky, QR, SVD

clear all; close all; clc

A = rand(100,100);
cond(A);
b = rand(100,1);

x = A\b;   % Matlab built in command (include sparse structure)
max(max(A*x-b));
norm(A*x-b);

A = rand(5,5);
[L,U] = lu(A);
norm(L*U-A);

[L,U,P] = lu(A);
norm(L*U-P*A);

A = hilb(7)   % doi xung + xac dinh duong
cond(A);
U = chol(A)   
norm(U' * U-A);

%% Example 5.9 - The role of condition number
err = [];  % vector sai so tuong doi

for n = 5:2:20
    n
    An = hilb(n);
    cond(An)
    
    xn = ones(n,1) ;   % nghiem chinh xac
    bn = An * xn;
    [L,U,P] = lu(An);
    y = L\(P*bn);
    x = U\y;
    err = [err norm(xn-x)/norm(xn)];
end

semilogy(5:2:20,err)
xlabel('n')
ylabel('Relative error')
legend('rel. err.')
grid on

%% QR and SVD

m = 10; n = 8;   % over-determined
A = rand(m,n)
[Q,R] = qr(A)


[U,S,V] = svd(A);

b = rand(m,1);
y = S\(U'*b);
x1 = V * y;
x2 = A\b
norm(x1-x2)


Aplus = pinv(A);

%% Image reduction using svd








