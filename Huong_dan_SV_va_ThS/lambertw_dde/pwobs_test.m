% Written by Shiming Duan, July. 31, 2010
% Dept. Mechanical Engineering, University of Michigan  
% Email: duansm@umich.edu                               
% For Matlab Toolbox of Time-Delay System Supplements   
% Piecewise observability test for a system of DDEs
% x_dot(t) = A*x(t) + Ad*x(t-T)
% y = C*x(t)
% Check the linear independence of the columns of C(sI-A-Ad*exp(-s*h))
% Let [O_1(s) O_2(s) ... O_n(s)]= C(sI-A-Ad*exp(-s*h)); If there exist p1, p2,...,pn, such that
% [O_1(p1) O_1(p1) ... O_n(p1);...;O_1(pn) O_2(pn) ...O_n(pn)] has full
% rank, then O_1(s), O_2(s), ..., O_n(s) are linearly independent

function  obs_flag = pwobs_test(A,Ad,C,h)

thred = 1e-5;
if nargin < 4
    error('No enough input arguments')
end

dim = length(A); % dimension of the system
syms s
char = s*eye(dim)-A-Ad*exp(-s*h); % char. eqn.
obs_gram = C*char^-1; 
runs = 10; % test 10 times. If failed everytime, the system is NOT piecewise observable; if passed for a run, the system is piecewise observable
obs_flag = 0; % default is set to be not piecewise observable

for j = 1:runs
    p = 10*rand(dim,1); % random variable [p1 p2 p3 ... pn]
    obs_matrix = [];
    for i = 1:dim
        obs_matrix = [obs_matrix;subs(obs_gram, p(i))]; % % formulating the matrix [O(p1);O(p2); .... ;O(pn)]
    end
    
    if det(obs_matrix) > thred % Code cua Duan co van de nang o day
        obs_flag =1;
        disp('The system of ODEs is piecewise observable.')
        return;
    end 
end

if obs_flag==0
    disp('Warning: The system of ODEs is NOT piecewise observable.')
end
    
    


