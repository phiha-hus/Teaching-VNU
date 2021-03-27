clear all; close all; clc;

nA = 4; nB=3; sP = nA+nB % total number of particles
t=2; k1=0.2; k2=0.2;

% rnl = randomNumberList
%rnl = [0.800, 0.801, 0.752, 0.661, 0.169, 0.956, 0.949,...
%    0.003, 0.201, 0.291, 0.615, 0.131, 0.241, 0.685, 0.116, 0.241, 0.849];

%n = length(rnl)
k1t = k1*t;
k2t=k2*t;

S = []; % Storing the number of nB after each time step

for j=1:2000 % Number of iteration until reaching t_max = 2000
for i=1: sP
    rnl = rand(1,1)
    if (rnl < nA/(nA+nB)) && (rnl<k1t) % This particle is A
        nA = nA-1;
        nB=nB+1;
    elseif (rnl >= nA/(nA+nB)) && (rnl < k2t) %This particle is B
        nA = nA+1;
        nB = nB-1;
    end
end

S = [S nB];
end

disp('The number of B-particle after t_max time step is')
nB

plot(S,'-')
%plot(linspace(1,2000,2000),S,'*')