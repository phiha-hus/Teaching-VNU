function Derec_AM2_DDAE11
%%-----------PP Adams Moulton 2 buoc, Bac 3 ap dung cho BT da bien doi VD11.=====================
clear all
format short e;
h0=input('nhap buoc thoi gian dau h0 = ');
%w=input('nhap w = ');
lamda = -1.5; w=10; a= 0.5;  b= 1; c = 0.8;
% ========== He so HEAB2, 2 buoc, order 2 ==================
beta0 = 5/12; beta1= 2/3; beta2= -1/12; 
% =============================== BDF 3 buoc, order 3 ====================
time = 20;
m=input('Nhap so lan giai m =');
tau = 1;
for i =1:m
 h(i)=2*h0/(2^i);
 % tong so diem chia la
 M = fix(time/h(i)); 
 % So diem chia trên mot khoang tre ====================
 no = fix(tau/h(i));
 % chia luoi thoi gia va gan gia tri chinh xac ============================
    for k = 1:M+no+3
        t(k)=(k- no -3)*h(i);
        VCX(k)= exp(lamda*t(k));   UCX(k)=(1 + w*t(k))*VCX(k);
    end;
 % gan gia tri ban dau ===================================================   
    for k=1:no+3
        U(k)= UCX(k);    V(k)=VCX(k); 
        SSU(k)= 0; SSV(k)= 0; 
    end;
 %========= GT ban dau cua Dao Ham EX====================
    DU(no+2)= (w + lamda + lamda*w*t(no+2))*exp(lamda*t(no+2)); DU(no+3)= (w + lamda + lamda*w*t(no+3))*exp(lamda*t(no+3));
    DV(no+2)= lamda*exp(lamda*t(no+2)); DV(no+3)= lamda*exp(lamda*t(no+3));
 %  tinh nghiem xap xy tai tung ðiem mot.=========================================
    for j=no+4:M+no+3
     % --- Dung noi suy---------------------------------
        VP1= U(j-1) - w*t(j)*V(j-1) + (beta1*DU(j-1) + beta2*DU(j-2))*h(i) - w*t(j)*(beta1*DV(j-1) + beta2*DV(j-2))*h(i) + a*beta0*(noisuy1(5,tau,t(j),h(i),t,V) - exp(lamda*t(j) - lamda*tau))*h(i);
        VP2= b*noisuy1(5,tau,t(j),h(i),t,U) + (c - b*w*t(j) + b*w*tau)*noisuy1(5,tau,t(j),h(i),t,V) - (b + c)*exp(lamda*t(j) - lamda*tau);
        V(j)= (VP1 - VP2*(1 - lamda*beta0*h(i)))/(1 - beta0*(w + lamda)*h(i)) ; 
        U(j)= (1 + w*t(j))*V(j) + VP2;
        DU(j)=(U(j) - U(j-1))/(beta0*h(i)) - (beta1*DU(j-1) + beta2*DU(j-2))/beta0;
        DV(j)=(V(j) - V(j-1))/(beta0*h(i)) - (beta1*DV(j-1) + beta2*DV(j-2))/beta0;
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
%plot(t,U,'- M');
title(' Do hoa nghiem xap xy X_1 (+) va nghiem CX X_1 (o)');
subplot(2,2,4);
plot(t,VCX);hold on
%plot(t,V,'- M');
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