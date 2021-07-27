%% Example 1.5 in Van Anh's thesis

clear all; close all; clc

% System coefficients
h = 0.1;                    % time-delay
tau = 0.0254; kp = 0.4451; ki = 2.3046; K = 1.53; 

A = [ 0 1; 0 -1/tau]; Ad = [0 0; -K*ki/tau -K*kp/tau] ; 

syms t
%g = [1; 1]; x0 = [0; 0];
g = [0; 0]; x0 = [1; 1];

vidu22_fun = @(t,y,z) A*y + Ad*z ;
vidu22_hist = @(t)[1;1] ;

t0 = 0; tf = 2; hstep = 1e-2; T = t0:hstep:tf;  

% NGHIEM CHINH XAC su dung dde23 thoi
sol = dde23(vidu22_fun,[h],vidu22_hist,[t0, tf]);
% sol = dde23(@vidu23_fun,[h],g,[t0, tf]);
sol_exact = deval(sol,T);

figure(1); clf;
subplot(2,2,1); 
plot(T,sol_exact)
title('Solution by dde23')
legend('x1','x2')

n_begin = 3; n_end = 3;
err = []; S_k = []; C_N = []; C_I = [];

for n = n_begin:1:n_end
    DDE_Sol = find_Sk(A,Ad,h,-n:n);
    sol_approx = [];
    S_k = []; C_N = []; C_I = [];
    
    for k = 1:(2*n+1)
        % disp(strcat('The matrix S_k coresponds to the branch: ',num2str(k-n-1),' is: '))
        Sk = DDE_Sol{k}.Sk;
        S_k(:,:,k) = Sk; 
        
        % CI tuong ung phan thuan nhat nên phu thuoc dieu kien ban dau
        [CI,~] = find_CI(A,Ad,h,Sk,g,x0);
        C_I(:,k) = CI;
        
        if k ==1
            alpha = real(max(eig(Sk)))
        end
    end 
    
    for i=1:length(T)
        t = T(i);
        xt = [0 ; 0];
        
        for k = 1:(2*n+1)
           temp1 = expm(S_k(:,:,k) .* t) * C_I(:,k); 
           xt = xt + temp1;                     
        end
        
        sol_approx = [sol_approx xt];
    end
    
    % err = [err abs(sol_exact - sol_approx)];      
    
end

subplot(2,2,2); 
plot(T,sol_approx)
title('Solution by LambertW')
legend('x1','x2')

% Vidu23_plot

subplot(2,2,3)
K = sol_exact(1)
sol_exact = sol_exact/K ;
plot(T,log(norm(sol_exact))./T)

subplot(2,2,4)
K = sol_approx(1)
sol_approx = sol_approx/K ;
plot(T,log(norm(sol_approx))./T)
