clear all: close all; clc 

a = -10:0.1:10;
b = -10:0.1:10;
c = -10:0.1:10;
[a,b,c] = meshgrid(a,b,c);

I = (a>0) & (c.*a-b.^2 > 0) & (4*b-2*a < 0) & ( -2*b.*a-24*c.*a+25*b.^2 + a.^2< 0) & (-8*c.*a+4*b .* c + 9 *b.^2+ 4 * c.^2 < 0);

scatter3(a(I),b(I),c(I),2)

xlabel('a')
ylabel('b')
zlabel('c')

view(40,40)		% change to top view
axis on