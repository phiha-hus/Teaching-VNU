r = -5:0.01:5;	% plotting range from -5 to 5
[x, y] = meshgrid(r);	% Get 2-D mesh for x and y based on r

cond1 = y-x.^2 > 0;	% check conditions for these values
cond2 = 4*x-2 < 0;
cond3 = -2*x-24*y+25*x.^2 < -1;
cond4 = -8*y+4*x .* y + 9 *x.^2+ 4 * y.^2 < 0;

output = ones(length(r)); % Initialize to 1
output(~(cond1 & cond2 & cond3 & cond4 )) = 0; % Zero out coordinates not meeting conditions.
imshow(output, 'xdata', r, 'ydata', r); % Display
xlabel('X'); ylabel('Y');
axis on;

