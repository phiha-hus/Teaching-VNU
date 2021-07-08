% Written by Shiming Duan, July. 31, 2010
% Dept. Mechanical Engineering, University of Michigan  
% Email: duansm@umich.edu                               
% For Matlab Toolbox of Time-Delay System Supplements
% DDE_Sol = find_Sk(A,Ad,h,N,Q_initial)
% The argument N and Q_initial are optional
% Find Sk for each branch given A, Ad and h for a DDE system 
% x_dot(t) = A*x(t) + Ad*x(t-T)
% Q_initial is the inital condition for Qk
% N is the index of the branches of interest; e.g., if one wants to calculate the branch
% {-2, -1 0, 1, 2,3}, then N = -2:3


function DDE_Sol = find_Sk(A,Ad,h,N,Q_initial,thredhold)
options_def =optimset('Display','iter');

if nargin < 3
    error('No enough input arguments')
elseif nargin < 4
    N = 0;    % default is to solve for principal branch
    Q_initial = expm(-A*h);
    options = options_def;
elseif nargin < 5 
    Q_initial = expm(-A*h);
    options = options_def;
elseif nargin <6
    options = options_def;
else
    if thredhold~=0 % if the thredhold is specified
        options = options_def;
        options.TolFun = thredhold; 
    else
        options = options_def;
    end
end

dim_Ad = length(Ad); % system order
DDE_Sol = {};
if dim_Ad == 1 % scalar case
    if Ad==0
        warning('First order ODE; x = exp(A*t)*x_0')
    else
        for i = 1:length(N)
             DDE_Sol{i}.Sk = lambertw(N(i),Ad*h*exp(-A*h))/h+A;
             DDE_Sol{i}.branchID = N(i);
             error1 = norm(lambertw(N(i),Ad*h*exp(-A*h))-Ad*h-h*DDE_Sol{i}.Sk);
             error2 = norm(DDE_Sol{i}.Sk*eye(length(A))-A-Ad*expm(-DDE_Sol{i}.Sk*h));
             DDE_Sol{i}.error = [error1 error2];
        end
    end
else % matrix case
       for i = 1:length(N)
             [DDE_Sol{i}.Sk,DDE_Sol{i}.Q,error1,error2] =QS_matrix(A,Ad,h,N(i),Q_initial,options);
             DDE_Sol{i}.branchID = N(i); 
             DDE_Sol{i}.error = [error1 error2];
       end
end

% currentQ =DDE_Sol{i}.Q % display the Q for current iteration



% find Sk and Qk for matrix case
function [S,Q1,error1,error2] = QS_matrix(A,Ad,h,branch_id,Q_initial,options)

if sum(sum(Q_initial == inf))>0
    reduced = 1; % Reduced elements for Q
    notfixed = Q_initial ~= inf; % find elements not fixed
    notfixed = reshape(notfixed,1,[]);
    res_Q = reshape(Q_initial,1,[]);
    reduced_Q = res_Q(notfixed); % formulate the new vecter of design variables
    Q = fsolve(@LambQ,reduced_Q,options,A,Ad,h,branch_id,reduced,notfixed);
else
    reduced = 0; % full elements for Q
    fixed = 0;
    Q = fsolve(@LambQ,Q_initial,options,A,Ad,h,branch_id,reduced,fixed);
end

if reduced ==0
    Q1 =Q;
else
    % if the reduced Q is used, convert the results back to matrix with
    % original dimension
    Q1 = zeros(1,length(A)^2); 
    position = find(notfixed==1);
    for k=1:length(position)
        Q1(position(k)) = Q(k);
    end
    Q1  = reshape(Q1,length(A),[]);
end      

H = Ad*h*Q1;
W = lambertw_matrix(branch_id,H);
S = W/h+A;
error1 = norm(W*expm(W+A*h)-Ad*h); % must be the solution to the lambert w function eqn.
error2 = norm(S*eye(length(A))-A-Ad*expm(-S*h)); % must be the solution to the char. eqn.
error = [error1 error2]; % both must be close to zero



