function NCE2HEMID_DDAE14
%-----------PP NCEHERK4 New: ERK order 4 + NCE order 2 ap dung cho BT da bien doi.=====================
% Công thuc NCE  duoc hieu chinh tu giai PT phi tuyen
% Tinh toan tren luoi deu h la uoc tau, tau = k.h.
clear all
format short e;
h0=input('nhap buoc thoi gian dau la uoc cua Pi, h0 = ');
time = 10*pi;
m=input('Nhap so lan giai m =');
tau = pi; 
thtakt=input('Vi tri nac kiem tra thuoc[0, 1], Thtakt = ');
% ========== He so c?a PP Trung diem hien, 2 buoc, order 2 =================
a21= 0.5;
b1= 0; b2= 1;
%==== He so cua noi suy ======================
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
 % chia luoi thoi gia va gan gia tri chinh xac ============================
 for k = 1:M+no+3
        t(k)=(k-no -3)*h(i);
        UCX(k)=exp(-t(k));   VCX(k)=sin(t(k));
        TT(k,1)=t(k); TT(k,2)= t(k) + h(i)/2; 
        % ---- Vi tri kiem tra va gia tri chinh xac tai diem KT
        TTKT(k)= t(k) + thtakt*h(i);
        UCXKT(k)=exp(-TTKT(k));   VCXKT(k)=sin(TTKT(k));
 end;
  % gan gia tri ban dau tu [- tau; 0] ===================================================   
 for k=1:no+2
    % Gia tri ban dau va sai so tai cac diem luoi
     U(k)= UCX(k); V(k)=VCX(k); SSU(k)= 0; SSV(k)= 0; 
    % Gia tri ban dau va sai so tai cac diem KT bac lien tuc 
     IVKT(k)= VCXKT(k); IUKT(k) = UCXKT(k); ErUC(k) = 0; ErVC(k)= 0; 
    % Gia tri noi suy tai cac diem nac ban dau
    IU(k,1)=exp(-TT(k,1));   IV(k,1)=sin(TT(k,1));
    IU(k,2)=exp(-TT(k,2));   IV(k,2)=sin(TT(k,2));
  end;
   % Them gia tri ban dau và sai so tai t = 0
   U(no+3)= UCX(no+3); V(no+3)=VCX(no+3); SSU(no+3)= 0; SSV(no+3)= 0; 
   %  tinh nghiem xap xy tai tung ðiem mot.=========================================
 for j=no+3:M+no+2
     %TT1 = t(j-1); TT2 = t(j-1)+ h(i)/2; TT3 = TT2; TT4 = t(j); 
     UTG1=U(j); VTG1=V(j);
     % Giai buoc thu nhat
     VP11=  UTG1*UTG1 + UTG1*VTG1*(TT(j,1)*TT(j,1) + 2*sin(TT(j,1))) + a21*UTG1*VTG1*(exp(-TT(j,1)) + 2*TT(j,1) + 2*cos(TT(j,1)))*h(i) + a21*UTG1*(sin(2*TT(j,1)))*h(i) + a21*(exp(-TT(j,1)))*((exp(-TT(j,1)))*IV(j-no,1) - exp(-TT(j,1)) + TT(j,1)*TT(j,1)*cos(TT(j,1)))*h(i);
     VP12= UTG1*(IV(j-no,2) + 1)*exp(-TT(j,2));
     VTG2 = (VP11 - VP12)/(UTG1*(TT(j,2)*TT(j,2) + 2*sin(TT(j,2)) + exp(-TT(j,2))));
     UTG2 = (VTG2 + IV(j-no,2) + 1)*exp(-TT(j,2));
     D1 =(UTG2 + VTG2*(TT(j,2)*TT(j,2) + 2*sin(TT(j,2))) - UTG1 - VTG1*(TT(j,1)*TT(j,1) + 2*sin(TT(j,1))))/(a21*h(i)); 
        % Giai buoc cuoi
     VP21= UTG2*UTG1 + UTG2*VTG1*(TT(j,1)*TT(j,1) + 2*sin(TT(j,1))) + UTG2*(b1*D1)*h(i) + b2*UTG2*VTG2*(exp(-TT(j,2)) + 2*TT(j,2) + 2*cos(TT(j,2)))*h(i) + b2*UTG2*(sin(2*TT(j,2)))*h(i) + b2*(exp(-TT(j,2)))*((exp(-TT(j,2)))*IV(j-no,2) - exp(-TT(j,2)) + TT(j,2)*TT(j,2)*cos(TT(j,2)))*h(i);
     VP22= UTG2*(IV(j-no+1,1) + 1)*exp(-t(j+1));
     V(j+1) = (VP21 - VP22)/(UTG2*(t(j+1)*t(j+1) + 2*sin(t(j+1)) + exp(-t(j+1))));
     U(j+1) = (V(j+1) + IV(j-no+1,1) + 1)*exp(-t(j+1));
     D2 =(U(j+1) + V(j+1)*(t(j+1)*t(j+1) + 2*sin(t(j+1))) - UTG1 - VTG1*(TT(j,1)*TT(j,1) + 2*sin(TT(j,1))))/(b2*h(i)) - (b1*D1)/b2;
         % Hieu chinh gia tri tai cac diem nac
     IV(j,1)=  (U(j) + V(j)*(TT(j,1)*TT(j,1) + 2*sin(TT(j,1))) + (b11*D1 + b21*D2)*h(i) - exp(-TT(j,1))*(IV(j-no,1) + 1))/(TT(j,1)*TT(j,1) + 2*sin(TT(j,1)) + exp(-TT(j,1)));
     IU(j,1)= exp(-TT(j,1))*(IV(j,1) +  IV(j-no,1) + 1);
     IV(j,2)=  (U(j) + V(j)*(TT(j,1)*TT(j,1) + 2*sin(TT(j,1))) + (b12*D1 + b22*D2)*h(i) - exp(-TT(j,2))*(IV(j-no,2) + 1))/(TT(j,2)*TT(j,2) + 2*sin(TT(j,2)) + exp(-TT(j,2)));
     IU(j,2)= exp(-TT(j,2))*(IV(j,2) +  IV(j-no,2) + 1);
     %IV(j,3)=  (U(j) + V(j)*(TT(j,1)*TT(j,1) + 2*sin(TT(j,1))) + (b13*D1 + b23*D2 + b33*D3 + b43*D4)*h(i) - exp(-TT(j,3))*(IV(j-no,3) + 1))/(TT(j,3)*TT(j,3) + 2*sin(TT(j,3)) + exp(-TT(j,3)));
     %%IU(j,3)= exp(-TT(j,3))*(IV(j,3) +  IV(j-no,3) + 1);
     %IV(j,4)=  (U(j) + V(j)*(TT(j,1)*TT(j,1) + 2*sin(TT(j,1))) + (b14*D1 + b24*D2 + b34*D3 + b44*D4)*h(i) - exp(-TT(j,4))*(IV(j-no,4) + 1))/(TT(j,4)*TT(j,4) + 2*sin(TT(j,4)) + exp(-TT(j,4)));
     %IU(j,4)= exp(-TT(j,4))*(IV(j,4) +  IV(j-no,4) + 1);
    %== gia tri tai diem nac KT bac lien tuc==================
     IVKT(j)=  (U(j) + V(j)*(TT(j,1)*TT(j,1) + 2*sin(TT(j,1))) + (bkt1*D1 + bkt2*D2)*h(i)  - exp(-TTKT(j))*(IVKT(j-no) + 1))/(TTKT(j)*TTKT(j) + 2*sin(TTKT(j)) + exp(-TTKT(j)));
     IUKT(j)=  exp(-TTKT(j))*(IVKT(j) +  IVKT(j-no) + 1);
     ErUC(j) = abs(UCXKT(j) - IUKT(j)); ErVC(j)= abs(VCXKT(j) - IVKT(j)); 
     % ---- Sai so ========================
     SSU(j+1)=abs(UCX(j+1)- U(j+1)); 
     SSV(j+1)=abs(VCX(j+1)- V(j+1));
 end;
  %Tinh sai so tai diem kiem tra deu
  EUktc(i)= max(ErUC);
  EVktc(i)= max(ErVC);
  %Tinh sai so bang Max tren toan mien
  cu(i) = max(SSU);   cv(i) = max(SSV);
  odh(i)=log(h(i)); ordu(i)= log(cu(i));  ordv(i)= log(cv(i));
 end
for e=1:(m-1)
    %Tinh bac sai so tai diem KT deu ---------
    orderCU(e)=log(EUktc(e)/EUktc(e+1))/log(2);
    orderCV(e)=log(EVktc(e)/EVktc(e+1))/log(2);
    %Tinh bac sai so bang Max tren toan mien
    orderU(e)=log(cu(e)/cu(e+1))/log(2);
    orderV(e)=log(cv(e)/cv(e+1))/log(2);
end;
%Tinh sai so bang Max tren toan mien va cap hoi tu roi rac 
disp('actual errors of U:'); disp(cu');
disp('Discrete orders of U:'); disp(orderU');
disp('actual errors of V:'); disp(cv');
disp('Discrete orders of V:'); disp(orderV');
%Tinh sai so tai diem giua va cap hoi tu deu
disp('Errors of U:'); disp(EUktc');
disp('Uniform order of U:'); disp(orderCU');
disp('Errors of V:'); disp(EVktc');
disp('Uniform order of V:'); disp(orderCV');
%=========== Ve do thi cap hoi tu roi rac cua U va V
subplot(2,2,1);
plot(odh,ordu, '+ M');hold off;
%plot(t,V);
title(' Do hoa cap hoi tu U');
subplot(2,2,2);
%plot(t,UCX);
plot(odh,ordv);
title(' Do hoa cap hoi tu V');
%subplot(2,2,3);
%plot(t,SSV);
%title(' Do hoa  sai so X_2 (-)');
%subplot(2,2,4);
%plot(t,SSU,'o M');
%title(' Do hoa sai so X_1 (o)');
% Lap hàm dau ra: Continuous output order 2 ================
%function kq = ConO2(tm,delay,stepSZ,DH,x)
%nu = floor((tm - delay)/stepSZ) + fix(delay/stepSZ) + 3; theta = (tm - delay)/stepSZ + 3 + fix(delay/stepSZ) - nu; 
%kq = x(nu) + stepSZ*(((2*theta)/3 - (theta^2)/2)*DH(nu,1) + (theta*DH(nu,2))/3 + (theta*DH(nu,3))/3 + ((theta^2)/2 - theta/3)*DH(nu,4));
 
% Lap hàm dau ra: Continuous output order 3 ================
%function kq = ConO3(tm,delay,stepSZ,DH,x)
%nu = floor((tm - delay)/stepSZ) + fix(delay/stepSZ) + 3; theta = (tm - delay)/stepSZ + 3 + fix(delay/stepSZ) - nu; 
%kq = x(nu) + stepSZ*((2*(theta^3)/3 - 3*(theta^2)/2 + theta)*DH(nu,1) + ((theta^2) - 2*(theta^3)/3)*DH(nu,2) + ((theta^2) - 2*(theta^3)/3)*DH(nu,3) + (2*(theta^3)/3 - (theta^2)/2 )*DH(nu,4));
