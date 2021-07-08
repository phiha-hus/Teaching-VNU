%% Vi du 2.2 phan Supplement

clear all; close all; clc

h     =  1;                                  % time-delay
A     =  [-1 -3;2 -5]; 
% Ad = A ;
Ad    =  [1.66 -0.697;0.93 -0.330]; 
B = [1 ; 0] ;

syms t
% g = [sin(t);cos(t)]; x0=[1;0];
% g = [1*t^0; 1*t^0];
g = [1; 0]; x0 = [1; 0];
% u = @(t) [cos(t); sin(t)];
u = @(t)cos(t) ;

vidu22_fun = @(t,y,z) A*y + Ad*z ;
vidu22_hist = @(t)[1;0] ;

t0 = 0; tf = 5; N = 1000; hstep = (tf-t0)/N; T = hstep .* [1:N];

% NGHIEM CHINH XAC su dung dde23 thoi
sol = dde23(vidu22_fun,[h],vidu22_hist,[t0, tf]);
% sol = dde23(@vidu23_fun,[h],g,[t0, tf]);
sol_exact = deval(sol,T);


%Q_ini = [inf inf;1 1];

% Neu chi so la -3:2 nen nhanh chinh k=0 thi phai tuong ung voi {4}
%Sk = DDE_Sol{4}.Sk ;

n_begin = 3;
n_end = 3;

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
        [CN,~] = find_CN(A,Ad,h,Sk);  
        C_N(:,:,k) = CN;
        
        % CI tuong ung phan thuan nhat nên phu thuoc dieu kien ban dau
        [CI,~] = find_CI(A,Ad,h,Sk,g,x0);
        C_I(:,k) = CI;
    end    
    
    for i=1:N
        t = T(i);
        xt = [0 ; 0];           
        
        for k = 1:(2*n+1)
           temp1 = expm(S_k(:,:,k) .* t) * C_I(:,k);
           
           % Dang bi loi o day 
           % Da sua lai bang cach viet ham tich phan moi quad_P
           temp2 = expm(S_k(:,:,k) .* t)  ;          
%            temp4 = C_N(:,:,k) * B ;
%            %temp3 = quad_P(@(s)expm(-s .* S_k(:,:,k)) * temp4 * u(s),0,t) ;        
%            syms s
%            temp3 = int(expm(-s .* S_k(:,:,k)) * C_N(:,:,k) * B * u(s),0,t) ;
           
           size(xt)
           size(temp1)
           size(temp2)
           size(temp3)
           
           xt = xt + temp1 + sum(temp2 * temp3);                     
         end
        
        sol_approx = [sol_approx xt];
    end
    
    err = [err abs(sol_exact - sol_approx)];      
    
    toc        
    
    % saiso = [saiso norm(sol_approx-sol_exact,Inf)]
end

Vidu23_plot