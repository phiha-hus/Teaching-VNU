function testGauss

n=3;

A = rand(n,n);
b = rand(n,1);

x = solve_Gauss(A,b)


% [L2,R2,P] = lu(A)
% err1 = norm(L2-inv(P) * L)
% err2 = norm(R2-R)
