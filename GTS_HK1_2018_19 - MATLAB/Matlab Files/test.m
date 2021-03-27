
n = 1e+3
S = 0;

tic
for i=1:n
  a = rand;
  S = S + a;
end
toc

X = [];
tic
for i=1:n  
  a = rand;
  X = [X a];
end
sum(X)
toc


