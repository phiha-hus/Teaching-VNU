clear all; close all; clc

v = -5:0.005:5;	% plotting range from -10 to 10
[x y] = meshgrid(v);	% get 2-D mesh for x and y

cond1 = y-x.^2 > 0;	% check conditions for these values
cond2 = 4*x-2 < 0;
cond3 = -2*x-24*y+25*x.^2 < -1;
cond4 = -8*y+4*x .* y + 9 *x.^2+ 4 * y.^2 < 0;

cond1 = double(cond1);	% convert to double for plotting
cond2 = double(cond2);
cond3 = double(cond3);
cond4 = double(cond4);

cond1(cond1 == 0) = NaN;	% set the 0s to NaN so they are not plotted
cond2(cond2 == 0) = NaN;
cond3(cond3 == 0) = NaN;
cond4(cond4 == 0) = NaN;

cond = cond1.*cond2.*cond3.*cond4;	% multiply the two condaces to keep only the common points
surf(x,y,cond)
colormap winter
grid on
xlim([-5,5]);
ylim([-5,5]);
view(0,90)		% change to top view
axis on


