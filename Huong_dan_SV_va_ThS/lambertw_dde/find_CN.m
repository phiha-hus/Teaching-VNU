% Written by Shiming Duan, Aug. 11, 2010
% Dept. Mechanical Engineering, University of Michigan  
% Email: duansm@umich.edu                               
% For Matlab Toolbox of Time-Delay System Supplements
% [CN, uniqueness]= find_CN(S,A,Ad,h)
% Find CN_k for the kth branch given A, Ad, h and Sk for this branch
% Important assumption: Any two eigenvalues from different branches are
% distinct !
% When Ad is invertible, S should be a matrix; when Ad is singular, S
% should be a scalar
% If the solution is unique, uniqueness=1; otherwise uniqueness=0

function [CN, uniqueness]= find_CN(A,Ad,h,S)
E = eig(S);
sys_order = length(A);

% scalar case
if sys_order ==1
    CN = 1/(1+Ad*h*exp(-S*h));
    uniqueness = 1;
    return
end


% matrix case
if length(S)==sys_order % general case, S is a matrix
    % S is a matrix here!
    syms s
    L1_num = 1;
    for i=1:sys_order
        L1_num = L1_num*(s-E(i));
    end
    char_eqn = det(s*eye(sys_order)-A-Ad*expm(-s*h)); % char eqn.
    L1 = diff(L1_num)/diff(char_eqn);
    R_matrix = [];
    L_matrix = [];
    for k =1:sys_order
        R1_num = adjugate(eye(sys_order)*E(k) - S);    
        L1_num = subs(L1,s,E(k));
        L2_num = adjugate(E(k)*eye(sys_order)-A-Ad*exp(-E(k)*h));
        R_matrix = [R_matrix;R1_num]; % stacked right side
        L_matrix = [L_matrix;L1_num*L2_num]; % stacked left side
    end
    uniqueness = rank(L_matrix)==sys_order; % for having an unique solution, matrix should have full column rank
    CN = pinv(R_matrix)*L_matrix; % use generalized pseudo inverse to calcuate the solution
elseif sys_order>1 && length(S)==1 % hybrid branch case, S is a scalar
    % S is a scalar here!
    syms s
    char_eqn = det(s*eye(sys_order)-A-Ad*expm(-s*h));
    L1 = 1/diff(char_eqn);
    L1_num = subs(L1,s,S);
    L2_num = adjugate(S*eye(sys_order)-A-Ad*exp(-S*h));
    CN = L1_num*L2_num;
    uniqueness = (L1_num~=0);
else
    error('Dimension of A and S are not compatible');
end
