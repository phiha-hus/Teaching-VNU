syms x
fplot(lambertw(x))
hold on
fplot(lambertw(-1,x))
hold off
axis([-0.5 4 -4 2])
legend('k=0','k=-1','Location','best')