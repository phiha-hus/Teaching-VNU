function HELM_DDAE14
%%-----------PP HE Adams Bashforth 2 buoc, Bac 2 hoac HELM3 ap dung cho BT da bien doi c?a VD14.=====================
clear all
format short e;
h0=input('nhap buoc thoi gian dau là uoc cua Pi, h0 = ');
% ========== He so HEAB2, 2 buoc, order 2 ==================
%beta1= 3/2; beta2= -1/2; 
% =============================== ELM3 buoc, order 3 ====================
beta1= 1/2; beta2= 3/2; beta3 = -1; 
time = 10*pi;%input('thoi gian la so nguyen < 5, time = ');
m=input('Nhap so lan giai m =');
tau = pi;
for i =1:m
 h(i)=2*h0/(2^i);
 % tong so diem chia la
 M = fix(time/h(i)); 
 % So diem chia trên mot khoang tre ====================
 no = fix(tau/h(i));
 % chia luoi thoi gia va gan gia tri chinh xac ============================
    for k = 1:M+no+3
        t(k)=(k-no -3)*h(i);
        UCX(k)= exp(-t(k));   VCX(k)= sin(t(k));
    end;
 % gan gia tri ban dau ===================================================   
    for k=1:no+3
        U(k)= UCX(k);    V(k)=VCX(k); 
        SSU(k)= 0; SSV(k)= 0; 
    end;
 %========= GT ban dau cua Dao Ham EX====================
    D(no+1)= -exp(-t(no+1)) + 2*t(no +1)*sin(t(no+1)) + t(no+1)*t(no +1)*cos(t(no+1)) + 2*sin(2*t(no+1)); 
    D(no+2)= -exp(-t(no+2)) + 2*t(no +2)*sin(t(no+2)) + t(no+2)*t(no +2)*cos(t(no+2)) + 2*sin(2*t(no+2));
 %  tinh nghiem xap xy tai tung ðiem mot.=========================================
    for j=no+4:M+no+3
        % --- Dung noi suy HEAB2 method---------------------------------
        %VP1= U(j-1)*U(j-1) + U(j-1)*V(j-1)*(t(j-1)*t(j-1) + 2*sin(t(j-1))) + h(i)*U(j-1)*(beta2*D(j-2) + beta1*sin(2*t(j-1))) + beta1*U(j-1)*V(j-1)*(exp(-t(j-1)) + 2*t(j-1) + 2*cos(t(j-1)))*h(i) + beta1*exp(-t(j-1))*(exp(-t(j-1))*(noisuy1(4,tau,t(j-1),h(i),t,V)) + t(j-1)*t(j-1)*cos(t(j-1)) - exp(-t(j-1)))*h(i);
        % --- Dung noi suy HELM3 method---------------------------------
        VP1= U(j-1)*U(j-1) + U(j-1)*V(j-1)*(t(j-1)*t(j-1) + 2*sin(t(j-1))) + h(i)*U(j-1)*(beta2*D(j-2) + beta3*D(j-3) + beta1*sin(2*t(j-1))) + beta1*U(j-1)*V(j-1)*(exp(-t(j-1)) + 2*t(j-1) + 2*cos(t(j-1)))*h(i) + beta1*exp(-t(j-1))*(exp(-t(j-1))*(noisuy1(4,tau,t(j-1),h(i),t,V)) + t(j-1)*t(j-1)*cos(t(j-1)) - exp(-t(j-1)))*h(i);
        VP2= U(j-1)*(noisuy1(4,tau,t(j),h(i),t,V) + 1)*exp(-t(j));
        V(j)= (VP1 - VP2)/(U(j-1)*(t(j)*t(j) + 2*sin(t(j)) + exp(-t(j))));
        U(j)= (1 + V(j) + noisuy1(4,tau,t(j),h(i),t,V))*exp(-t(j));
         % --- Dung noi suy HEAB2 method---------------------------------
        %D(j-1)=(U(j) + V(j)*(t(j)*t(j) + 2*sin(t(j))) - U(j-1) - V(j-1)*(t(j-1)*t(j-1) + 2*sin(t(j-1))))/(beta1*h(i)) - (beta2*D(j-2))/beta1;
         % --- Dung noi suy HELM3 method---------------------------------
        D(j-1)=(U(j) + V(j)*(t(j)*t(j) + 2*sin(t(j))) - U(j-1) - V(j-1)*(t(j-1)*t(j-1) + 2*sin(t(j-1))))/(beta1*h(i)) - (beta2*D(j-2) + beta3*D(j-3))/beta1;
        SSU(j)=abs(UCX(j)-U(j));
        SSV(j)=abs(VCX(j)-V(j));
    end;
    %Tinh sai so t?i diem cuoi
    %as(i) = SSU(M+no+2);
    %cs(i) = SSV(M+no+2);
    %Tinh sai so bang Max tren toan mien
    cu(i) = max(SSU);
    cv(i) = max(SSV);
end
for e=1:(m-1)
    %Tinh sai so t?i diem cuoi
    %orderU(e)=log(as(e)/as(e+1))/log(2);
    %orderV(e)=log(cs(e)/cs(e+1))/log(2);
    %Tinh sai so bang Max tren toan mien
    orderU(e)=log(cu(e)/cu(e+1))/log(2);
    orderV(e)=log(cv(e)/cv(e+1))/log(2);    
end;
disp('actual errors of U:'); disp(cu');
disp('error orders of U:'); disp(orderU');
disp('actual errors of V:'); disp(cv');
disp('error orders of V:'); disp(orderV');
subplot(2,2,1);
plot(t,SSU);
%plot(t,U,'- M');
title(' Do hoa sai so X_1');
subplot(2,2,2);
plot(t,SSV);
%plot(t,V,'- M');
title(' Do hoa sai so X_2');
subplot(2,2,3);
plot(t,UCX);hold on
plot(t,U,'- M');
title(' Do hoa nghiem xap xy X_1 (+) va nghiem CX X_1 (o)');
subplot(2,2,4);
plot(t,VCX);hold on
plot(t,V,'- M');
title(' Do hoa nghiem xap xy X_2 (+) va nghiem CX X_2 (o)');
% Lap cac ham NOI SUY
function KQ=noisuy1(sbns,kns,tns,h,t,y)
% sbns = so diem noi suy; kns = tau, la khoang noi suy; noisuy1(k,tau,Tns,h,t,U) = U(Tns - tau).
%m1=fix(kns/h); mht = fix((tns-t(1))/h)+ 1; m0 = mht-m1; 
m0 = fix((tns - kns - t(1))/h);
t0 = (tns- kns - t(m0))/h;
 for i=1:sbns-1
     D(i)= y(m0+i)- y(m0 + i -1);
 end
 ytg = y(m0)+ t0*D(1);
 tmoi=t0;
 for j =2:sbns-1
     for i =1: sbns-j
         D(i)= D(i+1)- D(i);
     end
     tmoi=tmoi*(t0 - j + 1)/j;
     ytg = ytg + tmoi*D(1);
 end
% D
 KQ=ytg;