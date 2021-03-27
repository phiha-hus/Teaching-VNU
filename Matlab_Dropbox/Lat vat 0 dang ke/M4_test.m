clear all; close all; clc

U1 = rand(1e3);

% a = 1; b = 10;
% U2 = (b-a) * U1 + a
% 
% 
% figure(1); clf;
% 
% subplot(1,2,1);
% hist(U2,25)
% title('Histogram')
% 
% subplot(1,2,2);
% U3 = -log(U1)/2;
% hist(U3,20)


t = 0.0; count = 0; 
dt = -log(U1)/2; 
lambda = 2;

while (max(t+dt)<1)
    count = count + 1;
    U1 = rand(1e3);
    dt = -log(U1)/2;
    t = t + dt
end





