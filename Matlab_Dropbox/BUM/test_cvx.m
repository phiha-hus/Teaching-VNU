N = 4096;
SRF = 16;
fc = N/(2*SRF);
s = 64;

Fn = dftmtx(N);
Fn1 = Fn(1:(fc+1), :);
Fn2 = Fn((N-fc+1):N, :);
Fn = [Fn2;Fn1];

dist = 31;
%
support = find_support(N,SRF,dist,s);
    
%u = rand(1,length(support));
x = zeros(N,1);
% x(support) = (-1).^randint(length(support),1) .* (10.^rand(length(support),1));
% x(support) = exp(1i* rand(1,length(support)));
[U,S,V] = svds(Fn(:, support));
[min_value, min_ind] = min(svds(Fn(:,support)));
x(support) = V(:,min_ind);

cvx_begin
    variable z(N) complex;
    Fn*z == Fn*x;
    minimize norm(z,1)
cvx_end
