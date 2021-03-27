clear all; close all; clc

v = -10:0.1:10;
[x,y] = meshgrid(v);  % create a grid
ineq1 = -x.^2 + y > 0;    % some inequality
f1 = double(ineq1);
surf(x,y,f1)
hold on
%view(0,90)            % rotate surface plot to top view

ineq2 = -4*x + 2 > 0;    % some inequality
f2 = double(ineq2);
surf(x,y,f2)
hold on

