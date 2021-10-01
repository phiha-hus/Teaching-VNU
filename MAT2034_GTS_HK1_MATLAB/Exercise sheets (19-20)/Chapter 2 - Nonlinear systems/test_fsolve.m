clear all; close all; clc

F = @(x) [2*x(1) - x(2) - exp(-x(1));
         -x(1) + 2*x(2) - exp(-x(2))];
     
x0 = [-5;-5];     

options = optimoptions('fsolve','Display','iter');

tic
[x,fval] = fsolve(F,x0,options);
toc


options = optimset('Display','iter');
options.MaxIter = 1e+3 ;
options.MaxFunEvals = 3e+3 ;
options.FunctionTolerance = 1e-12;
options.OptimalityTolerance = 1e-12;
options.StepTolerance = 1e-12;

tic
[x,fval] = fsolve(F,x0,options);
toc