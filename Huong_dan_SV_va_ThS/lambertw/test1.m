
x = linspace(-6,6,51);
y = linspace(-6,6,51);
[x,y] = meshgrid(x,y);
w = lambertw(x + i*y);

surf(x,y,abs(w));
axis([-6,6,-6,6,0,2.5]);
view(40,32);
xlabel('Re z');
ylabel('Im z');
title('|W_0(z)|');