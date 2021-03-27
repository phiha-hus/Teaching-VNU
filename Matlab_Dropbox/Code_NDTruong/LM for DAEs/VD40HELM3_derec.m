function VD40HELM3_derec
clear all
format short e;
% === Cac PP HELM3, HEAB3, HEAB2 truc tiep ================================
beta1 = 0.5; beta2 = 1.5; beta3 = -1; % HELM3; p = 2
%beta1 = 23/12; beta2 = -16/12; beta3 = 5/12; % Adams Bashforth 3; p =3
%beta1 = 3/2; beta2 = -1/2; % Adams Bashforth 2; p =2
h0 = input('nhap buoc thoi gian dau h0 = ');
m=7; %input('Nhap so lan giai m =');
time = 100;
for k =1:m
    h(k)=2*h0/(2^k);
    N = time/h(k);
    for i=1:N+1
        t(i)=(i-1)*h(k);  UCX(i)=exp(-t(i)); VCX(i)= sin(t(i));
    end;
    %U xap xi cua x1, V xap xi cua x2 (neu 2 buoc thi bo gia tri xuat phat thu 3)
    U(1)=1; V(1)=0; U(2)= UCX(2); V(2)= VCX(2); U(3)= UCX(3); V(3)= VCX(3); 
    SSU(1)=0; SSV(1)=0; SSU(2)= 0; SSV(2)= 0; SSU(3)= 0; SSV(3)= 0;
    DU(1)= -1; DV(1)= 1; DU(2)= - exp(-t(2)); DV(2)= cos(t(2)); 
    % Neu 2 buoc thi for i = 2: N ============================
    for i=3:N
        % HELM 3 step
        VP =  U(i)*U(i) + U(i)*V(i)*(t(i)*t(i)+ 2*sin(t(i)) + h(k)*beta1*exp(-t(i))) + U(i)*(beta2*DU(i-1) + beta3*DU(i-2))*h(k) + U(i)*(t(i)*t(i)+ 2*sin(t(i)))*(beta2*DV(i-1) + beta3*DV(i-2))*h(k) + h(k)*beta1*(U(i)*sin(2*t(i)) - V(i)*exp(-2*t(i))) + h(k)*beta1*exp(-t(i))*(t(i)*t(i)*cos(t(i)) - exp(-t(i)));
        % ============== Adams Bashforth 2; p =2
         %VP =  U(i)*U(i) + U(i)*V(i)*(t(i)*t(i)+ 2*sin(t(i)) + h(k)*beta1*exp(-t(i))) + U(i)*(beta2*DU(i-1))*h(k) + U(i)*(t(i)*t(i)+ 2*sin(t(i)))*(beta2*DV(i-1))*h(k) + h(k)*beta1*(U(i)*sin(2*t(i)) - V(i)*exp(-2*t(i))) + h(k)*beta1*exp(-t(i))*(t(i)*t(i)*cos(t(i)) - exp(-t(i)));
         %-----------------------------
        V(i+1) = (U(i)*(sin(t(i+1)) - 1)*exp(-t(i+1)) + VP)/(U(i)*(exp(-t(i+1)) + t(i)*t(i) + 2*sin(t(i))));
        U(i+1) = V(i+1)*exp(-t(i+1)) - (exp(-t(i+1)))*(sin(t(i+1)) - 1); 
        % --------- HELM 3 step
        DU(i)= (U(i+1)  - U(i))/(beta1*h(k)) - (beta2*DU(i-1) + beta3*DU(i-2))/beta1; 
        DV(i)= (V(i+1)  - V(i))/(beta1*h(k)) - (beta2*DV(i-1) + beta3*DV(i-2))/beta1;
        % Adams Bashforth 2; p =2
        %DU(i)= (U(i+1)  - U(i))/(beta1*h(k)) - (beta2*DU(i-1))/beta1; 
        %DV(i)= (V(i+1)  - V(i))/(beta1*h(k)) - (beta2*DV(i-1))/beta1;
        SSU(i+1)=abs(U(i+1)-UCX(i+1));    SSV(i+1)=abs(V(i+1)-VCX(i+1));
    end;
    a(k)=max(SSU);
    b(k)=max(SSV);
end
for e=1:(m-1)
    orderCXU(e)=log(a(e)/a(e+1))/log(2);
    orderCXV(e)=log(b(e)/b(e+1))/log(2);
end;
%c=a';
%d=order';
%D=[b; c]
subplot(2,2,1);
plot(t(1:(N+1)),U(1:(N+1)),'-')
title ('nghiem xap xy x1 (-), Ham e')
subplot(2,2,2);
plot(t(1:(N+1)),V(1:(N+1)),'-')
title ('nghiem xap xy x2 (-), Ham sin')
subplot(2,2,3);
plot(t(1:(N+1)),UCX(1:(N+1)),'-')
title ('nghiem dung (-), Ham e')
subplot(2,2,4);
plot(t(1:(N+1)),VCX(1:(N+1)),'-')
title ('nghiem dung (-), Ham sin')
disp('actual errors of U:'); disp(a');
disp('error orders of U:'); disp(orderCXU');
disp('actual errors of V:'); disp(b');
disp('error orders of V:'); disp(orderCXV');