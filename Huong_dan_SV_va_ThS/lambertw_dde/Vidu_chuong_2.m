%% Vi du 2.3
clear all; close all; clc

h     =  1;                                  % time-delay
A     =  [-1 -3;2 -5]; 
Ad    =  [1.66 -0.697;0.93 -0.330]; 

syms t
% g = [sin(t);cos(t)]; x0=[1;0];
% g = [1*t^0; 1*t^0];
g = [0; 0]; x0 = [1; 1];
u = @(t) [0; 0];

t0 = 0; tf = 5; N = 1000; hstep = (tf-t0)/N; T = hstep .* [1:N];

% NGHIEM CHINH XAC su dung dde23 thoi
sol = dde23(@vidu23_fun,[h],@vidu23_hist,[t0, tf]);
% sol = dde23(@vidu23_fun,[h],g,[t0, tf]);
sol_exact = deval(sol,T);
plot(T,sol_exact)

%Q_ini = [inf inf;1 1];

% Neu chi so la -3:2 nen nhanh chinh k=0 thi phai tuong ung voi {4}
%Sk = DDE_Sol{4}.Sk ;

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
        % C_N(:,k) = CN;
        
        % CI tuong ung phan thuan nhat nên phu thuoc dieu kien ban dau
        [CI,flag] = find_CI(A,Ad,h,Sk,g,x0);
        C_I(:,k) = CI;
    end
    
    % Something wrong with the code
    % S_k computing using DDE_Sol is somehow different 
    S_k(:,:,1) = [-0.3499 + 4.9801i  -1.6253 + 0.1459i; 2.4174 + 0.1308i  -5.1048 + 4.5592i];
    S_k(:,:,3) = [-0.3499 - 4.9801i  -1.6253 - 0.1459i; 2.4174 - 0.1308i  -5.1048 - 4.5592i];
   
    
    for i=1:N
        t = T(i);
        xt = 0;   
        
        for k = 1:(2*n+1)
            temp1 = expm(S_k(:,:,k) .* t) * C_I(:,k); 
            % xt = temp1 + sum(temp2 .* integral(@(s)exp(-s .* S_k) .* u(s),0,t,'ArrayValued',true));         
            xt = xt + temp1;
        end
        
        sol_approx = [sol_approx xt];
    end
    
    err = [err abs(sol_exact - sol_approx)];      
    
    toc        
    
    % saiso = [saiso norm(sol_approx-sol_exact,Inf)]
end


Vidu23_plot