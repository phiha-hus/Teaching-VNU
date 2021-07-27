%% Vi du 2.1   
clear all; close all; clc

a = -1; a_d = 0.5; h = 1;
temp = a_d * h * exp(-a * h); 

syms t
g = 1 * t^0; x0 = 1; 

n = 15;

clf; figure(1);
grid on

for k = -n:1:n
    k ;
    s_k  = 1/h * lambertw(k,temp) + a ;
    plot(s_k,'r *'); hold on;
    grid on
    xlim([-6,2])
    xlabel('Re(s_k)')
    ylabel('Im(s_k)')
    legend('s_k')   
end

%print('hinh2.1','-dpng')

K = []; S_k = []; C_i = []; C_n = [];

for k = 0:1:n
    K = [K k] ;
    
    s_k  = 1/h * lambertw(k,temp) + a ; 
    S_k = [S_k s_k];
    
    den = 1 + a_d * h * exp(-s_k * h);
        
    % D_ik = ( x0 + a_d * exp(-s_k * h ) * integral( @(t)exp(-s_k * t ) * 1,0,h) )/den ;
    [C_ik, uniqueness] = find_CI(a,a_d,h,s_k,g,x0);
    
    % Ki?m tra xem code l?p trình có trùng v?i k?t qu? trong tool không.
    % D_ik - C_ik
    
    C_i = [C_i C_ik];
    %C_i = [C_i [D_ik;C_ik]];
    
    C_nk = 1/den; 
    C_n = [C_n C_nk];
end

% [K; S_k ; C_i; C_n]

S_k 

n = 5;
DDE_Sol = find_Sk(a,a_d,h,-n:n);
for i = 1:11
    tic
    disp(strcat('The matrix S_k coresponds to the branch: ',num2str(i-6),' is: '))
    DDE_Sol{i}.Sk
    toc
end

%=========================================================================
%% Vi du 2.2 

A = -1; Ad = 0.5; h = 1; B = 1; 

syms t
g = [1]; x0 = [1];

u = @(t)sin(t);

t0 = 0; tf = 5; N = 1000; 
hstep = (tf-t0)/N; T = hstep .* [1:N];

n_begin = 3; n_end = 3;
err = [];

%=========================================================================
% Tim nghiem dua vao dde23
vidu22_hist = @(t) 1;
vidu22_fun = @(t,y,z) A * y + Ad * z + sin(t);
sol = dde23(vidu22_fun,h,vidu22_hist,[t0,tf]);
sol_exact = deval(sol,T);
%=========================================================================

for n = n_begin:2:n_end

    tic
    
    DDE_Sol = find_Sk(A,Ad,h,-n:n);
    
    sol_approx = []; 
    
    S_k = []; C_N = []; C_I = [];
    
    for k = 1:(2*n+1)
        % disp(strcat('The matrix S_k coresponds to the branch: ',num2str(k-n-1),' is: '))
        Sk = DDE_Sol{k}.Sk;        
        S_k = [S_k Sk];

        % CI tuong ung phan thuan nhat nên phu thuoc dieu kien ban dau
        % Minh tinh ra CI khac han trong bai bao la sao nhi 
        % [CI,flag] = find_CI(A,Ad,h,Sk,g,x0);
        
        temp1 = 1 + Ad * h * exp(-Sk * h);
        temp2 = integral(@(t) exp(-Sk * (t+h)),0,h);
        CI = (x0 + Ad * temp2)/temp1 ;
        C_I = [C_I CI];
        
        % CN tuong ung phan u(t) nen 0 phu thuoc dieu kien ban dau
        [CN,flag] = find_CN(A,Ad,h,Sk)  ;        
        C_N = [C_N CN];        
    end
    
    C_I
    C_I = [0.0016+0.0005i  0.0038+0.0015i 0.0197+0.0111i 0.9422 0.0197-0.0111i 0.0038-0.0015i 0.0016-0.0005i]
    
  
    for i=1:N
        t = T(i);       
        temp1 = exp(S_k * t) * C_I';                
        temp2 = exp(S_k * t) .* C_N;        
        xt = temp1 + sum(temp2 .* integral(@(s)exp(-s .* S_k) .* u(s),0,t,'ArrayValued',true));                         
        sol_approx = [sol_approx xt];
    end     
    
     saiso = abs(sol_exact - sol_approx);    
     err = [err norm(saiso,Inf)];      

     figure(n); clf;
     subplot(2,3,1)
     plot(T,sol_approx,'r -.','LineWidth',2)
     legend('Sol. approx.')
     grid on
     
     subplot(2,3,2)
     semilogy(T,saiso,'r -.','LineWidth',2);
     ylim([-0.5 13])
     legend('Abs. err.')
     grid on    
     
     subplot(2,3,3)
     semilogy(T,saiso./abs(sol_approx),'b.-','LineWidth',2); 
     ylim([-0.5 13])
     legend('Rel. err.')
     grid on 
     
     print(gcf, '-dpdf', 'Vidu_2_2.pdf');
     
    toc    
end

% clf; figure(10);
% semilogy(n_begin:1:n_end,err(1,:),'r *'); 
% grid on
% xlabel('So thanh phan xap xi n')
% ylabel('Sai so tuyet doi')
% legend('1')   
%=========================================================================

