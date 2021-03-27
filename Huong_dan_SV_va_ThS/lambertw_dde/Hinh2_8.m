% Plot figure 2.8
clear all; close all; clc

clf; figure(1);

S = [-0.12 + 181.88i, -0.12 - 181.88i, -176.73 - 428.66i, -176.73 + 428.66i, -91.61, -233.30 - 755.05i, -233.30 + 755.05i];

real_S = real(S);
imag_S = imag(S);

plot(real_S,imag_S,'r *')
xlim([-250 20])
ylim([-1000 1000])
xlabel('Re')
ylabel('Im')
grid on