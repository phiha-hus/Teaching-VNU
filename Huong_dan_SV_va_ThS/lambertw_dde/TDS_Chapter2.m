% Written by Sun Yi Sept. 30, 2009
% Dept. Mechanical Engineering, University of Michigan  
% Phone: (734) 763-2227
% Email: syjo@umich.edu                               
% Matrix Function: Compute Qk

%%
global h k1 k2 A Ad                 % Define variable for fsolve
 
k1 = 1; k2=k1;                          % Branch of the matrix Lambert W function
I     =  eye(2);                          % for identity matrix
h     =  1;                                  % time-delay
A     =  [-1 -3;2 -5]; Ad    =  [1.66 -0.697;0.93 -0.330]; % coefficient matrices of the system A=-Ad;B=-A
x0=I;                                             % initial condition for interation => identity matrix  

tol=1e-8; options = optimset('MaxFunEvals',1000,'TolFun',1e-8); % iteration options
X   =  fsolve(@equations,x0,options)                          % iteration

%% 

Q     =  [X(1) X(2);X(3) X(4)]                                  % result for Qk
D     =  Ad * h * Q;[v,d] =  eig(D);                                % To compute Sk    
W     =  v*[lambertw(k1,d(1,1)) 0;0 lambertw(k2,d(2,2))]*inv(v);
Sk    =  1/h * W + A                                                % result for Sk  

%% Phan nay em khong can quan tam, truoc het cu lam phan tren lien quan den Sk di cai da  
% 

[V_Sk, D_Sk] = eig(Sk);                                 % eigenvalue of Sk

eig_Sk = [D_Sk(1,1) D_Sk(2,2)];

    for jj=1:2

        s=eig_Sk(jj);

        %%%%%%%%%%%%%%%%%%second

        up=2*s-eig_Sk(1)-eig_Sk(2);                                               % Equations (2.58) and (2.60) in book 

        down=2*s+307/125*exp(-s)+133/100*s*exp(-s)+6-10041/50000*exp(-s)^2;

        Second=up/down;

        %%%%%%%%%%%%%%%%%%third

        Char=s*I-Ad*exp(-s*T)-A;

        Third=[Char(2,2) -Char(1,2);-Char(2,1) Char(1,1)];

        %%%%%%%%%%%%%%%%%%Fourth

        Fourth=[1;0]+Ad*[1/s-exp(-s*h)/s;0];                                  % Using Initial Conditions

        %%%%%%%%%%%%%%%%%%left

        Left(:,jj)=Second*Third*Fourth;

        Left1(:,:,jj)=Second*Third;

        %%%%%%%%%%%%%%%%%%first

        First(:,:,jj)=s*I+[-Sk(2,2) Sk(1,2);Sk(2,1) -Sk(1,1)];

    end

    %%%%%%%%%%%%%%%%%% Ck

    for ii=1:1

        Rig=[First(ii,1,1) First(ii,2,1);First(ii,1,2) First(ii,2,2)];

        C_I = inv(Rig)*[Left(ii,1);Left(ii,2)];

    end

    %%%%%%%%%%%%%%%%%% Ck prime

    Rig1=[First(1,1,1) First(1,2,1);First(1,1,2) First(1,2,2)];

    Cprime = inv(Rig1)*[Left1(1,1,1) Left1(1,2,1);Left1(1,1,2) Left1(1,2,2)];

