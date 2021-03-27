% Matrix lambert W function
W = [1 1 1;0 1 2;0 0 3]; branch = -1; % repeated eigenvalue
jordan(W)
lambertw_matrix(branch,W)

W1 = [1 1 1;0 1 2;0 0 0]; branch = 3; % hybrid branch
jordan(W1)
lambertw_matrix(branch,W1)


%%
% find Sk, CI and CN
% #1 
h     =  1;                                  % time-delay
A     =  [-1 -3;2 -5]; Ad    =  [1.66 -0.697;0.93 -0.330]; % coefficient matrices of the system A=-Ad;B=-A
Sk = [-1.9676-23.5358i -1.6379+0.0245i;2.1921+0.0593i  -6.6492-23.4771i];
[CN,flag] = find_CN(A,Ad,h,Sk)
syms t
g = [sin(t);cos(t)];x0=[1;0];
g = [1*t^0; 1*t^0];
[CI,flag] = find_CI(A,Ad,h,Sk,g,x0)


%%
% find Sk, CI and CN
% #2
h     =  1;                                  % time-delay
A = [0 1;-9.397/2 0];
Ad = [0 0 ;0 0.5];
Q_ini = [inf inf;1 1];
DDE_Sol = find_Sk(A,Ad,h,-3:2,Q_ini);
DDE_Sol{1}.Sk
Sk = eig(DDE_Sol{1}.Sk)
[CN,flag] = find_CN(A,Ad,h,Sk(1));
syms t
g = [sin(t);0];x0=[1;0];
[CI,flag] = find_CI(A,Ad,h,Sk(1),g,x0)

%%
% find Sk, CI and CN
%#3
A     =  [-1 0;1 -2]; Ad    =  [0.1 0;0 0]; h=1;
Q_ini = [1 1;inf inf];
DDE_Sol = find_Sk(A,Ad,h,-2,Q_ini);
DDE_Sol{1}.Sk
Sk = eig(DDE_Sol{1}.Sk)
[CN,flag] = find_CN(A,Ad,h,Sk(1));
syms t
g = [sin(t-h);cos(t)];x0=[1;0];
[CI,flag] = find_CI(A,Ad,h,Sk(1),g,x0);

%%
% find Sk, CI and CN
% #4

h = 0.91;m=0;T0=408000;N=480;np=0.70;c=4.3;delta=1.57;k=c/(N*T0);  % Parameters
A=  [-delta 0 0;(1-np)*N*delta -c 0;np*N*delta 0 -c]                                       % Coefficients in state-space forms
Ad= [0 k*T0*exp(-m*h) 0;0 0 0;0 0 0]
Q_ini   =  [inf inf inf;1.0e+003*(-2.8456-7.3652i)  1.0e+003*(0.0218+0.1126i) 0;inf inf inf];
DDE_Sol = find_Sk(A,Ad,h,1,Q_ini,1e-15); % rightmost eigenvalue is contained in branch 1. You may check branch -1, 0, 1 to find this rightmost eigenvalue. 
Sk = eig(DDE_Sol{1}.Sk)
[CN,flag] = find_CN(A,Ad,h,Sk(1));
syms t
g = [sin(t);cos(t);1];x0=[1;0;0];
[CI,flag] = find_CI(A,Ad,h,Sk(1),g,x0);

%%
% find Sk, CI and CN
% #5
A = -0.1; Ad = -0.1; h=1; 
DDE_Sol = find_Sk(A,Ad,h,-2);
Sk = DDE_Sol{1}.Sk;
[CN,uniqueness] = find_CN(A,Ad,h,Sk)
syms t
g = [sin(t)];x0=[5];
[CI,uniqueness] = find_CI(A,Ad,h,Sk,g,x0)


%%

% contr. and obs. test
h     =  1;                                  % time-delay
A = [0 1;-9.397 0];
Ad = [0 0 ;-2 -3];
C = [0 1];
pwobs_test(A,Ad,C,h)

A = [-1 0;1 -2]; Ad = [0.1 0.0001;0 0]; h = 1; B = [0;1];
pwcontr_test(A,Ad,B,h)
%%

% pole placement
h     =  0.2;                                  % time-delay
A = [0 1;-9.397/2 0]; Ad = [0 0 ;0 0];B = [0;1]; % open loop
% pole_desired=[-1+2.5i];
pole_desired=[-1+2i];
% pole_desired=[-1+1.5i];
% Q_ini = [1.0709-0.0000i 0.2916-0.0000i;-0.2310-0.0000i 0.9697-0.0000i];
% Q_ini = [inf inf;-0.2310-0.0000i 0.9697-0.0000i];
Q_ini = [inf inf;1 1];
Kd_ini = [1 1];
contr_mode = 2; % u = Kd*x(t-h)
branch =0;
% K  = place_dde(A,Ad,B,h,branch,pole_desired,contr_mode);
threshold = [1e-10, 1e-8];
K = place_dde(A,Ad,B,h,branch,pole_desired,contr_mode,Kd_ini,Q_ini,threshold);


%%
% observer design
h = 0.91;m=0;T0=408000;N=480;np=0.70;c=4.3;delta=1.57;k=c/(N*T0);  % Parameters
A=  [-delta 0 0;(1-np)*N*delta -c 0;np*N*delta 0 -c]                                       % Coefficients in state-space forms
Ad= [0 k*T0*exp(-m*h) 0;0 0 0;0 0 0]
C = [0 1 1];
Q_ini   =  [inf inf inf;1.0e+003*(-2.8456-7.3652i)  1.0e+003*(0.0218+0.1126i) 0;inf inf inf];

L_ini = [0 0 0]';
pole_desired = -1.2;
contr_mode = 4; % observer mode
branch = 1;
K = place_dde(A,Ad,C,h,branch,pole_desired,contr_mode,L_ini,Q_ini);

%%

% stability radiu #1
I=eye(2);B=[0;1];h=0.1;
K=[-0.6971   -1.6893];Kd=[-0.6971   -1.6893];% -0.5 -6
A=[0 0;0 1]+B*K; Ad=[-1 -1;0 -0.9]+B*Kd;
sr = stabilityradius_dde(A,Ad,I,I,I,h) % E1=I, E2=I, F=I
%%
% stability radiu #2
I=eye(2);B=[0;1];h=0.1;
A=[0 0;0 1]; Ad=[-1 -1;0 -0.9];
sr = stabilityradius_dde(A,Ad,I,I,I,h) % E=I, F1=I, F2=I
%%

% gramian test
load gramian_test.mat
t1 = 4; % observability at t1 = 4 sec.
C = [0 1];
B = [0;1];
ob_gramm_lambert = obs_gramian_dde(Result(1:4),C,t1)
det(ob_gramm_lambert)
ct_gramm_lambert = contr_gramian_dde(Result(1:4),B,t1)

