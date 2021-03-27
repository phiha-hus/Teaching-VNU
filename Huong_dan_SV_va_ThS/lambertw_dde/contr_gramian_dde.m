% Written by Shiming Duan, Aug. 06, 2010
% Dept. Mechanical Engineering, University of Michigan  
% Email: duansm@umich.edu                               
% For Matlab Toolbox of Time-Delay System Supplements   
% Calculate the controllability gramian at time t1 given the solution set of a DDE
% Sol.Sk{i} and Sol.CN{i} are the S and CN for the ith branch
% It is important that the the two conjugated branches must be used or not used the same time 
% t1 is the time instant of the gramian
% B is the matrix for input u(t) = B*u(t)
function gramian = contr_gramian_dde(Sol,B,t1)
syms xi
total_branch = length(Sol);
kernal_lambert = {};
gramm_lambert = {};
det_lambert =[];
sum_e = zeros(length(Sol{1}.CN),1);

% approximate the kernal func using the first n branches
for j=1:total_branch
    sum_e = sum_e+expm(Sol{j}.Sk*xi)* Sol{j}.CN*B;
end
kernal_lambert  = sum_e*sum_e';
% approximate the gramian func using the first n branches
gramian = vpa(int(kernal_lambert,0,t1));
