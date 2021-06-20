%%

clear all; close all; clc

%% contr. and obs. test
clear all; close all; clc
h =  1;                                  % time-delay
A = [-27 3.6 6 ; 9.6 -12.5 0 ; 0 9 -5];
Ad = [0 0 0 ; 21 0 0 ;0 0 0];
B = [0.26 0 ; -0.9 -0.8 ; 0 0.18];

% pwcontr_test(A,Ad,B,h) -> Fail for nonsquare contr_Gramian 

n = 1;
DDE_Sol = find_Sk(A,Ad,h,-n:1:n);

S_k = []; E_k = [];

for i = 1:1:(2*n+1)
    k = i - n - 1;
    %disp(strcat('The matrix S_k coresponds to the branch: ',num2str(k),' is: '))
    S_k = DDE_Sol{i}.Sk ;
    E_k = [E_k eig(S_k)];
end

E_k

clf; figure(1);
plot(E_k,'r *')
xlim([-100 10])
grid on

%%
clear all; close all; clc
h =  1;                                  % time-delay
A = [-27 3.6 6 ; 9.6 -12.5 0 ; 0 9 -5];
%Ad = [0 0 0 ; 21 0 0 ;0 0 0];
Ad = zeros(3,3) ;

B = [0.26 0 ; -0.9 -0.8 ; 0 0.18]; 
C = [0 1 0];

% B = [0.26; -0.9; 0]; % Code cua Duan co van de nang o day
% C = [1 0 0; 0 1 0];

disp('Using the test pwcontr_test: ')
contr_flag = pwcontr_test(A,Ad,B,h)

[~,n] = size(A);

disp('Ctrb of the pair (A,B): ')
contr_flag = (rank(ctrb(A,B))==n)

syms s w
% exp(-s*h)
dim = length(A); % dimension of the system
char = s*eye(dim)-A-Ad*w ;
det(char);
contr_gram =  (char^-1 .* det(char)) * B ;
simplify(contr_gram) ;  % Cha tac dung may

rank(contr_gram) ;

[~,nB] = size(B);
contr_flag = (rank(contr_gram) == nB);
if contr_flag ==1;
    disp('The system of ODEs is piecewise controllable.')
else
    disp('The system of ODEs is NOT piecewise controllable.')
end


pwobs_test(A,Ad,C,h)

%%


%% gramian test
% load gramian_test.mat
% t1 = 4; % observability at t1 = 4 sec.
% C = [0 1];
% B = [0;1];
% ob_gramm_lambert = obs_gramian_dde(Result(1:4),C,t1)
% det(ob_gramm_lambert)
% ct_gramm_lambert = contr_gramian_dde(Result(1:4),B,t1)
