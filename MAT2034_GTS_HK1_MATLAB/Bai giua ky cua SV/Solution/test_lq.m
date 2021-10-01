T = [1 2 3 4 5 6 7]';
V = [2.31 2.01 1.80 1.66 1.55 1.47 1.41]';

n = length(T);

A = [ones(n,1) T T.^2];

value = (A' * A ) \ (A'*V);
flipdim(value)'

disp('So sanh ket qua voi polyfit')
polyfit(T,V,2)