%% Vi du 2.3
clear all; close all; clc

h     =  1;                                  % time-delay
A     =  [-1 -3;2 -5]; 
Ad    =  [1.66 -0.697;0.93 -0.330]; 

syms t
% g = [sin(t);cos(t)]; x0=[1;0];

% Example in Paper 2003
g = [1; 0]; x0 = [1; 0];

% Example in Paper 2014
% g = [0; 0]; x0 = [1; 1];
u = @(t) [0; 0];

t0 = 0; tf = 5; N = 5000; hstep = (tf-t0)/N; T = hstep .* [1:N];

% NGHIEM CHINH XAC su dung dde23 thoi
sol = dde23(@vidu23_fun,[h],@vidu23_hist,[t0, tf]);
% sol = dde23(@vidu23_fun,[h],g,[t0, tf]);
sol_exact = deval(sol,T);


%%  Tim S_k trong ham Lambert

%Q_ini = [inf inf;1 1];

n_begin = 1;
n_end = 1;

err = [];
S_k = []; C_N = []; C_I = [];

for n = n_begin:1:n_end

    DDE_Sol = find_Sk(A,Ad,h,-n:n);
    
    sol_approx = [];
    
    tic    
    
    S_k = []; C_N = []; C_I = [];
    
    for k = 1:(2*n+1)
        % disp(strcat('The matrix S_k coresponds to the branch: ',num2str(k-n-1),' is: '))
        Sk = DDE_Sol{k}.Sk;
        S_k(:,:,k) = Sk;

        % CN tuong ung phan u(t) nen 0 phu thuoc dieu kien ban dau
        % [CN,flag] = find_CN(A,Ad,h,Sk)  
        % C_N(:,:,k) = CN;
        
        % CI tuong ung phan thuan nhat nên phu thuoc dieu kien ban dau
        [CI,flag] = find_CI(A,Ad,h,Sk,g,x0);
        C_I(:,k) = CI;
    end
    
%     S_k 
%     
%     % S_k computing using DDE_Sol is somehow different 
%     S_k(:,:,1) = [-0.3499 - 4.9801i  -1.6253 + 0.1459i ; 2.4174 + 0.1308i  -5.1048 - 4.5592i] 
%     S_k(:,:,2) = [0.3055 + 0.0000i  -1.4150 - 0.0000i  ; 2.1317 + 0.0000i  -3.3015 + 0.0000i] 
%     S_k(:,:,3) = [-0.3499 + 4.9801i  -1.6253 - 0.1459i ; 2.4174 - 0.1308i  -5.1048 + 4.5592i]    
   
    
    % Approximate using 7 terms
%    S_k(:,:,1) = [-1.6528-17.2617i  -1.6392-0.0337i ; 2.2017-0.0800i  -6.3426-17.1805i ];
%    S_k(:,:,2) = [-1.1893-11.0140i  -1.6416-0.0549i ; 2.2340-0.1191i  -5.9012-10.8791i];
%    S_k(:,:,3) = [-0.3499-4.9801i  -1.6253-0.1459i ; 2.4174-0.1308i  -5.1048-4.5592i];
%    S_k(:,:,4) = [0.3055+0.0000i  -1.4150-0.0000i ;  2.1317+0.0000i  -3.3015+0.0000i];
%    S_k(:,:,5) = [-0.3499+4.9801i  -1.6253+0.1459i ; 2.4174+0.1308i  -5.1048+4.5592i] ;
%    S_k(:,:,6) = [-1.1893+11.0140i  -1.6416+0.0549i ; 2.2340+0.1191i  -5.9012+10.8791i];
%    S_k(:,:,7) = [-1.6528+17.2617i  -1.6392+0.0337i ; 2.2017+0.0800i  -6.3426+17.1805i ] ;
     
      
   C_I
   % CI different between the code and the result presented in the paper
   % using Laplace transforms
   C_I = [1.3663+3.9491i -1.7327 1.3663-3.9491i ; 3.2931+9.3999i -6.5863 3.2931-9.3999i]
   
   
    for i=1:N
        t = T(i);
        xt = 0;   
        
        for k = 1:(2*n+1)
            temp1 = expm(S_k(:,:,k) .* t) * C_I(:,k);               
            xt = xt + temp1;
        end
        
        sol_approx = [sol_approx xt];
    end
    
    err = [err abs(sol_exact - sol_approx)];      
    
    toc        
    
    % saiso = [saiso norm(sol_approx-sol_exact,Inf)]
end

Vidu23_plot

%% Vi du 2.3 theo cach khac 
% Using only 3 terms and use already built matrices 

clear all; close all; clc

h     =  1;                                  % time-delay
A     =  [-1 -3;2 -5]; 
Ad    =  [1.66 -0.697;0.93 -0.330]; 

syms t
% g = [sin(t);cos(t)]; x0=[1;0];

% Example in Paper 2003
g = [1; 0]; x0 = [1; 0];

% Example in Paper 2014
% g = [0; 0]; x0 = [1; 1];

u = @(t) [0; 0];

t0 = 0; tf = 5; N = 1000; hstep = (tf-t0)/N; T = hstep .* [1:N];

% NGHIEM CHINH XAC su dung dde23 thoi
sol = dde23(@vidu23_fun,[h],@vidu23_hist,[t0, tf]);
%sol = dde23(@vidu23_fun,[h],g,[t0, tf]);
sol_exact = deval(sol,T);


S_k(:,:,1) = [-0.3499 + 4.9801i  -1.6253 + 0.1459i; 2.4174 + 0.1308i  -5.1048 + 4.5592i];
S_k(:,:,2) = [0.3055 + 0.0000i  -1.4150 + 0.0000i; 2.1317 + 0.0000i  -3.3015 + 0.0000i];
S_k(:,:,3) = [-0.3499 - 4.9801i  -1.6253 - 0.1459i; 2.4174 - 0.1308i  -5.1048 - 4.5592i];

C_I = [1.3663+3.9491i -1.7327 1.3663-3.9491i ; 3.2931+9.3999i -6.5863 3.2931-9.3999i];


x = @(t) expm(S_k(:,:,1) .* t) * C_I(:,1) + expm(S_k(:,:,2) .* t) * C_I(:,2) + expm(S_k(:,:,3) .* t) * C_I(:,3);

% Dont understant while error "Matrix dimensions must agree".
% Loi nay sinh o day 
% sol_approx = x(T);

for i  = 1:N
     t = T(i);
     sol_approx(:,i) =  x(t);
end


err = [];
err = [err abs(sol_exact - sol_approx)];      

Vidu23_plot