clear all; close all; clc

x = -10:0.1:10;
y = -10:0.1:10;
[x,y] = meshgrid(x,y);
z = y-x.^2;
surf(x,y,z)
colorbar