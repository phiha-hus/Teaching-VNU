function NCEHEMID_DDAE11
%-----------PP NCEHEMID New: ERK order 2 + NCE order 2 ap dung cho BT da bien doi.=====================
% NCE-2  duoc hieu chinh tu giai PT phi tuyen 
% Tinh toan tren luoi deu h là uoc cua tau, tau = no.h
clear all
format short e;
h0=input('nhap buoc thoi gian dau la uoc cua tau, h0 = ');
w=10; lamda= -1.5; a=0.5; b=1; c=0.8;
time = 20;
m=input('Nhap so lan giai m =');
tau = 1; 
thtakt=input('Vi tri nac kiem tra thuoc[0, 1], Thtakt = ');
% ========== He so c?a PP Explicit Midpoint, 2 buoc, order 2 ===========6=======
a21= 0.5; %a32= 0.5; a43 = 1; 
b1= 0; b2= 1; %b3= 1/3; b4 = 1/6;
theta1 = 0; theta2 = 0.5;
%== He so NCE2 uniform order 2 ===========
b11=  theta1 - (theta1^2); b12 = theta2 - (theta2^2);
b21= theta1^2; b22 = theta2^2;
%== Tinh he so noi suy tai diem KT bac deu NCE2 ===========
bkt1= thtakt - thtakt^2; bkt2= thtakt^2;
% ====== CHUONG TRINH TINH TOAN CHINH ====================
for i =1:m
 h(i)=2*h0/(2^i);
 % tong so diem chia la
 M = time/h(i); 
 % So diem chia trên mot khoang tre ====================
 no = tau/h(i);
 %==== He so cua noi suy ======================
 atg =  tau/(h(i))- no;
 % chia luoi thoi gia va gan gia tri chinh xac ============================
 for k = 1:M+no+3
        t(k)=(k-no -3)*h(i);
        VCX(k)=exp(lamda*t(k));   UCX(k)=(1 + w*t(k))*VCX(k);
        TT(k,1)=t(k); TT(k,2)=t(k) + h(i)/2; 
        % ---- Vi tri kiem tra va gia tri chinh xac tai diem KT
        TTKT(k)= t(k) + thtakt*h(i);
        VCXKT(k)=exp(lamda*TTKT(k));   UCXKT(k)=(1 + w*TTKT(k))*VCXKT(k);
 end;
  % gan gia tri ban dau tu [- tau; 0] ===================================================   
 for k=1:no+2
    % Gia tri ban dau va sai so tai cac diem luoi
     U(k)= UCX(k);    V(k)=VCX(k); SSU(k)= 0; SSV(k)= 0; 
     % Gia tri noi suy tai cac diem nac ban dau
    iT(k,1)= t(k) + theta1*h(i); iT(k,2)= iT(k,1) + h(i)/2; 
    IV(k,1)=exp(lamda*iT(k,1));   IU(k,1)=(1 + w*iT(k,1))*IV(k,1);
    IV(k,2)=exp(lamda*iT(k,2));   IU(k,2)=(1 + w*iT(k,2))*IV(k,2);
    IVKT(k)= VCXKT(k); IUKT(k) = UCXKT(k); ErUC(k) = 0; ErVC(k)= 0; 
  end;
   % Them gia tri ban dau và sai so tai t = 0
   U(no+3)= UCX(no+3); V(no+3)=VCX(no+3); SSU(no+3)= 0; SSV(no+3)= 0; 
   %  tinh nghiem xap xy tai tung ðiem mot.=========================================
 for j=no+3:M+no+2
     TT1 = t(j); TT2 = t(j)+ h(i)/2;  
     % Giai buoc thu nhat
     UTG1=U(j); VTG1=V(j);
     % Giai buoc thu hai
     VP11=  UTG1 - w*VTG1*TT1 + a21*lamda*(UTG1 - w*VTG1*TT1)*h(i) + a21*a*(IV(j-no,1) - exp(lamda*TT1 - tau*lamda))*h(i);
     VP12= b*IU(j-no,2) + (c - b*w*TT2 + b*w*tau)*IV(j-no,2) - (b + c)*exp(lamda*TT2 - lamda*tau);
     VTG2 = VP11 - VP12; UTG2 = VP12 + (w*TT2 + 1)*VTG2;
     D1 =(UTG2 - w*VTG2*TT2 - UTG1 + w*VTG1*TT1)/(a21*h(i)); 
      % Giai buoc cuoi
     VP21=  UTG1 - w*VTG1*TT1 + (b1*D1)*h(i) + b2*lamda*(UTG2 - w*VTG2*TT2)*h(i) + b2*a*(IV(j-no,2) - exp(lamda*TT2 - lamda*tau))*h(i);
     VP22= b*IU(j-no +1,1) + (c - b*w*t(j+1) + b*w*tau)*IV(j-no +1,1) - (b + c)*exp(lamda*t(j+1) - lamda*tau);
     V(j+1)= VP21 - VP22;   U(j+1)= VP22 + (w*t(j+1) + 1)*V(j+1); 
      D2 =(U(j+1) - w*V(j+1)*t(j+1) - UTG1 + w*VTG1*TT(j,1))/(b2*h(i));
     % Hieu chinh gia tri tai cac diem nac
     iT(j,1)= t(j) + theta1*h(i); iT(j,2)= iT(j,1) + h(i)/2; 
     IV(j,1)= U(j) - w*V(j)*t(j) + (b11*D1 + b21*D2)*h(i) - b*IU(j-no,1) - (c + b*w*tau - b*w*iT(j,1))*IV(j-no,1) + (b + c)*exp(lamda*iT(j,1) - tau*lamda);
     IU(j,1)= (1 + w*iT(j,1))*IV(j,1) +  b*IU(j-no,1) + (c + b*w*tau - b*w*iT(j,1))*IV(j-no,1) - (b + c)*exp(lamda*iT(j,1) - tau*lamda);
     IV(j,2)= U(j) - w*V(j)*t(j) + (b12*D1 + b22*D2)*h(i) - b*IU(j-no,2) - (c + b*w*tau - b*w*iT(j,2))*IV(j-no,2) + (b + c)*exp(lamda*iT(j,2) - tau*lamda);
     IU(j,2)= (1 + w*iT(j,2))*IV(j,2) +  b*IU(j-no,2) + (c + b*w*tau - b*w*iT(j,2))*IV(j-no,2) - (b + c)*exp(lamda*iT(j,2) - tau*lamda);
       %== gia tri tai diem nac KT bac lien tuc==================
     IVKT(j)=  U(j) - w*V(j)*t(j) + (bkt1*D1 + bkt2*D2)*h(i) - b*IUKT(j-no) - (c + b*w*tau - b*w*TTKT(j))*IVKT(j-no) + (b + c)*exp(lamda*TTKT(j) - tau*lamda);
     IUKT(j)= (1 + w*TTKT(j))*IVKT(j) +  b*IUKT(j-no) + (c + b*w*tau - b*w*TTKT(j))*IVKT(j-no) - (b + c)*exp(lamda*TTKT(j) - tau*lamda);
     ErUC(j) = abs(UCXKT(j) - IUKT(j)); ErVC(j)= abs(VCXKT(j) - IVKT(j)); 
     % ---- Sai so ========================
     SSU(j+1)=abs(UCX(j+1)- U(j+1));
     SSV(j+1)=abs(VCX(j+1)- V(j+1));
 end;
 %Tinh sai so tai diem kiem tra deu
  EUktc(i)= max(ErUC);
  EVktc(i)= max(ErVC);
  %Tinh sai so bang Max tren toan mien
  cu(i) = max(SSU);
  cv(i) = max(SSV);   
end
for e=1:(m-1)
    %Tinh bac sai so tai diem KT deu ---------
   orderCU(e)=log(EUktc(e)/EUktc(e+1))/log(2);
    orderCV(e)=log(EVktc(e)/EVktc(e+1))/log(2);
    %Tinh bac sai so bang Max tren toan mien
    orderU(e)=log(cu(e)/cu(e+1))/log(2);
    orderV(e)=log(cv(e)/cv(e+1))/log(2);
end;
%Tinh sai so tai diem giua
disp('Uniform Errors of U:'); disp(EUktc');
disp('Uniform order of U:'); disp(orderCU');
disp('Uniform Errors of V:'); disp(EVktc');
disp('Uniform order of V:'); disp(orderCV');
%Tinh sai so bang Max tren toan mien, 
disp('actual errors of U:'); disp(cu');
disp('Discrete orders of U:'); disp(orderU');
%Tinh sai so t?i diem cuoi
disp('actual errors of V:'); disp(cv');
disp('Discrete orders of V:'); disp(orderV');
subplot(2,2,1);
plot(t,VCX, '+ M');
plot(t,V);
title(' Do hoa nghiem xap xy X_2 (o)');
subplot(2,2,2);
plot(t,UCX);
plot(t,U);
title(' Do hoa nghiem xap xy X_1 ');
subplot(2,2,3);
plot(t,SSV);
title(' Do hoa  sai so X_2 (-)');
subplot(2,2,4);
plot(t,SSU,'o M');
title(' Do hoa sai so X_1 (o)');
% Lap hàm dau ra: Continuous output order 2 ================
%function kq = ConO2(tm,delay,stepSZ,DH,x)
%nu = floor((tm - delay)/stepSZ) + fix(delay/stepSZ) + 3; theta = (tm - delay)/stepSZ + 3 + fix(delay/stepSZ) - nu; 
%kq = x(nu) + stepSZ*(((2*theta)/3 - (theta^2)/2)*DH(nu,1) + (theta*DH(nu,2))/3 + (theta*DH(nu,3))/3 + ((theta^2)/2 - theta/3)*DH(nu,4));
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
