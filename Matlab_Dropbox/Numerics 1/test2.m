% Comment: Time difference while computing only svd and both svd and rank
% of the same matrix is REALLY NOT SMALL.

m = 1e2;

T1 = []; T2 =[];
for i=1:20
    n = m*i;
A = rand(n,n);
tic
 [U,S,V]=svd(A);
 ra = rank(S);
T1 = [T1 toc];

tic
 [U,S,V]=svd(A);
T2 = [T2 toc];
end

T1-T2