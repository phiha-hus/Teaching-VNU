% Written by Shiming Duan, July. 31, 2010
% Dept. Mechanical Engineering, University of Michigan  
% Email: duansm@umich.edu                               
% For Matlab Toolbox of Time-Delay System Supplements   
% Piecewise controllability test for a system of DDEs
% x_dot(t) = A*x(t) + Ad*x(t-T)+B*u(t)
% Check the linear independence of the rows of (sI-A-Ad*exp(-s*h))B
% Let [C_1(s);C_2(s); ... ;C_n(s)]= (sI-A-Ad*exp(-s*h))B; If there exist p1, p2,...,pn, such that
% [C_1(p1) C_1(p2) ... C_1(pn);...;C_n(p1) C_n(p2) ...C_n(pn)] has full
% rank, then C_1(s), C_2(s), ... ,C_n(s) are linearly independent

function contr_flag = pwcontr_test(A,Ad,B,h)

thred = 1e-5;
if nargin < 4
    error('No enough input arguments')
end

dim = length(A); % dimension of the system
syms s
char = s*eye(dim)-A-Ad*exp(-s*h);  % Day la ham so

contr_gram = char^-1*B;

[~,nB] = size(B);

runs = 10; % test 10 times. If failed everytime, the system is NOT piecewise controllable; if passed for a run, the system is piecewise controllable
contr_flag = 0; % default is set to be not piecewise controllable

for j = 1:runs
    p = 10*rand(dim,1); % random variable [p1 p2 p3 ... pn]
    
    contr_matrix = [];
    for i = 1:dim
        contr_matrix = [contr_matrix subs(contr_gram, p(i))]; % formulating the matrix [C(p1) C(p2) .... C(pn)]
    end
    
    rank(contr_matrix)
    
    % if det(contr_matrix) > thred    % Code cua Duan co van de nang o day
     if rank(contr_matrix) == dim    
        contr_flag =1;
        disp('The system of ODEs is piecewise controllable.')
        break;
    end 
end

if contr_flag==0
    disp('The system of ODEs is NOT piecewise controllable.')
end
    
    


