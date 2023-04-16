f1 = @(x,y) x.^2 .* sin(y);
disp(integral2(f1,0,4,0,pi));
f2 = @(x,y) x.^2 .*(x + y);
x1 = @(y) y;
disp(integral2(f2,0,2,x1,3));
