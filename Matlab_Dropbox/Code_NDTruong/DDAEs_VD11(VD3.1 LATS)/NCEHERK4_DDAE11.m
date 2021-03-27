function NCEHERK4_DDAE11
%-----------PP NCEHERK-4 New: ERK order 4 + NCE order 3 ap dung cho BT da bien doi.=====================
% Công thuc NCE  duoc hieu chinh tu giai PT phi tuyen
% Tinh toan tren luoi deu, buoc luoi h la uoc tau hay tau = no.h
clear all
format short e;
h0=input('nhap buoc thoi gian dau la uoc cua tau, h0 = ');
w=10; lamda= -1.5; a=0.5; b=1; c=0.8;
time = 20;
m=input('Nhap so lan giai m =');
tau = 1; 
thtakt=input('Vi tri nac kiem tra thuoc[0, 1], Thtakt = ');
% ========== He so c?a PP ERK, 4 buoc, order 4 ===========6=======
a21= 0.5; a32= 0.5; a43 = 1; b1= 1/6; b2= 1/3; b3= 1/3; b4 = 1/6;
%==== He so cua noi suy ======================
theta2 = 0.5;
%== NCE uniform order 2 ===========
b11= 0; b12= (2*theta2)/3 - (theta2^2)/2; b13 = b12; b14 = 2/3 - 1/2;
b21= 0; b22= theta2/3; b23 = b22; b24 = 1/3;
b31= 0; b32= theta2/3; b33 = b32; b34 = 1/3;
b41= 0; b42= (theta2^2)/2 - theta2/3; b43 = b42; b44 = 0.5 - 1/3;
%== Tinh he so noi suy tai diem KT bac deu NCE2 ===========
bkt1= (2*thtakt)/3 - (thtakt^2)/2; bkt2= thtakt/3; bkt3 = bkt2; bkt4 = (thtakt^2)/2 - thtakt/3;
%== NCE uniform order 3 ===========
%b11= 0; b12= (2*(theta2^3))/3 - 3*(theta2^2)/2 + theta2; b13 = b12; b14 = 2/3 - 3/2 + 1;
%b21= 0; b22= theta2^2 - 2*(theta2^3)/3; b23 = b22; b24 = 1 - 2/3;
%b31= 0; b32= theta2^2 - 2*(theta2^3)/3; b33 = b32; b34 = 1 - 2/3;
%b41= 0; b42= 2*(theta2^3)/3 - (theta2^2)/2; b43 = b42; b44 = 2/3 - 0.5;
%== Tinh he so noi suy tai diem KT bac deu NCE3 ===========
%bkt1= (2*(thtakt^3))/3 - 3*(thtakt^2)/2 + thtakt; bkt2= thtakt^2 - 2*(thtakt^3)/3; bkt3 = bkt2; bkt4 = 2*(thtakt^3)/3 - (thtakt^2)/2;
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
        VCX(k)=exp(lamda*t(k));   UCX(k)=(1 + w*t(k))*VCX(k);
        TT(k,1)=t(k); TT(k,2)=t(k) + h(i)/2; TT(k,3)= TT(k,2); TT(k,4)= t(k) + h(i);
        % ---- Vi tri kiem tra va gia tri chinh xac tai diem KT
        TTKT(k)= t(k) + thtakt*h(i);
        VCXKT(k)=exp(lamda*TTKT(k));   UCXKT(k)=(1 + w*TTKT(k))*VCXKT(k);
 end;
  % gan gia tri ban dau tu [- tau; 0] ===================================================   
 for k=1:no+2
    % Gia tri ban dau va sai so tai cac diem luoi
     U(k)= UCX(k);    V(k)=VCX(k); SSU(k)= 0; SSV(k)= 0; 
    % Gia tri ban dau va sai so tai cac diem KT bac lien tuc 
     IVKT(k)= VCXKT(k); IUKT(k) = UCXKT(k); ErUC(k) = 0; ErVC(k)= 0; 
    % Gia tri noi suy tai cac diem nac ban dau
    IV(k,1)=exp(lamda*TT(k,1));   IU(k,1)=(1 + w*TT(k,1))*IV(k,1);
    IV(k,2)=exp(lamda*TT(k,2));   IU(k,2)=(1 + w*TT(k,2))*IV(k,2);
    IV(k,3)=IV(k,2);   IU(k,3)=IU(k,2);
    IV(k,4)=exp(lamda*TT(k,4));   IU(k,4)=(1 + w*TT(k,4))*IV(k,4);
    % Gia tri dao ham tai cac diem nac ban dau
    %DH(k,1)= (lamda + w + lamda*w*TT(k,1))*exp(lamda*TT(k,1));
    %DH(k,2)= (lamda + w + lamda*w*(TT(k,2)))*exp(lamda*(TT(k,2))); DU(k,3) = DU(k,2); 
    %DH(k,4)= (lamda + w + lamda*w*(TT(k,4)))*exp(lamda*(TT(k,4)));
 end;
   % Them gia tri ban dau và sai so tai t = 0
   U(no+3)= UCX(no+3); V(no+3)=VCX(no+3); SSU(no+3)= 0; SSV(no+3)= 0; 
   %  tinh nghiem xap xy tai tung ðiem mot.=========================================
 for j=no+3:M+no+2
     %TT1 = t(j-1); TT2 = t(j-1)+ h(i)/2; TT3 = TT2; TT4 = t(j); 
     UTG1=U(j); VTG1=V(j);
     % Giai buoc thu nhat
     VP11=  UTG1 - w*VTG1*TT(j,1) + a21*lamda*(UTG1 - w*VTG1*TT(j,1))*h(i) + a21*a*(IV(j-no,1) - exp(lamda*TT(j,1) - tau*lamda))*h(i);
     VP12= b*IU(j-no,2) + (c - b*w*TT(j,2) + b*w*tau)*IV(j-no,2) - (b + c)*exp(lamda*TT(j,2) - lamda*tau);
     VTG2 = VP11 - VP12; UTG2 = VP12 + (w*TT(j,2) + 1)*VTG2;
     D1 =(UTG2 - w*VTG2*TT(j,2) - UTG1 + w*VTG1*TT(j,1))/(a21*h(i)); 
      % Giai buoc thu hai
     VP21=  UTG1 - w*VTG1*TT(j,1) + a32*lamda*(UTG2 - w*VTG2*TT(j,2))*h(i) + a32*a*(IV(j-no,2) - exp(lamda*TT(j,2) - lamda*tau))*h(i);
     VP22= b*IU(j-no,3) + (c - b*w*TT(j,3) + b*w*tau)*IV(j-no,3) - (b + c)*exp(lamda*TT(j,3) - lamda*tau);
     VTG3 = VP21 - VP22; UTG3 = VP22 + (w*TT(j,3) + 1)*VTG3;
     D2 =(UTG3 - w*VTG3*TT(j,3) - UTG1 + w*VTG1*TT(j,1))/(a32*h(i)); 
      % Giai buoc thu ba
     VP31=  UTG1 - w*VTG1*TT(j,1) + a43*lamda*(UTG3 - w*VTG3*TT(j,3))*h(i) + a43*a*(IV(j-no,3) - exp(lamda*TT(j,3) - lamda*tau))*h(i);
     VP32= b*IU(j-no,4) + (c - b*w*TT(j,4) + b*w*tau)*IV(j-no,4) - (b + c)*exp(lamda*TT(j,4) - lamda*tau);
     VTG4 = VP31 - VP32; UTG4 = VP32 + (w*TT(j,4) + 1)*VTG4;
     D3 =(UTG4 - w*VTG4*TT(j,4) - UTG1 + w*VTG1*TT(j,1))/(a43*h(i)); 
      % Giai buoc cuoi
     VP41=  UTG1 - w*VTG1*TT(j,1) + (b1*D1 + b2*D2 + b3*D3)*h(i) + b4*lamda*(UTG4 - w*VTG4*TT(j,4))*h(i) + b4*a*(IV(j-no,4) - exp(lamda*TT(j,4) - lamda*tau))*h(i);
     VP42= b*IU(j-no,4) + (c - b*w*t(j+1) + b*w*tau)*IV(j-no,4) - (b + c)*exp(lamda*t(j+1) - lamda*tau);
     V(j+1)= VP41 - VP42;   U(j+1)= VP42 + (w*t(j+1) + 1)*V(j+1); 
      D4 =(U(j+1) - w*V(j+1)*t(j+1) - UTG1 + w*VTG1*TT(j,1))/(b4*h(i)) - D1 - 2*D2 - 2*D3;
     % Hieu chinh gia tri tai cac diem nac
     IV(j,1)=  U(j) - w*V(j)*TT(j,1) - b*IU(j-no,1) - (c + b*w*tau - b*w*TT(j,1))*IV(j-no,1) + (b + c)*exp(lamda*TT(j,1) - tau*lamda);
     IU(j,1)= (1 + w*TT(j,1))*IV(j,1) +  b*IU(j-no,1) + (c + b*w*tau - b*w*TT(j,1))*IV(j-no,1) - (b + c)*exp(lamda*TT(j,1) - tau*lamda);
     IV(j,2)=  UTG1 - w*VTG1*TT(j,1) + (b12*D1 + b22*D2 + b32*D3 + b42*D4)*h(i) - b*IU(j-no,2) - (c + b*w*tau - b*w*TT(j,2))*IV(j-no,2) + (b + c)*exp(lamda*TT(j,2) - tau*lamda);
     IU(j,2)= (1 + w*TT(j,2))*IV(j,2) +  b*IU(j-no,2) + (c + b*w*tau - b*w*TT(j,2))*IV(j-no,2) - (b + c)*exp(lamda*TT(j,2) - tau*lamda);
     IV(j,3)=  UTG1 - w*VTG1*TT(j,1) + (b13*D1 + b23*D2 + b33*D3 + b43*D4)*h(i) - b*IU(j-no,3) - (c + b*w*tau - b*w*TT(j,3))*IV(j-no,3) + (b + c)*exp(lamda*TT(j,3) - tau*lamda);
     IU(j,3)= (1 + w*TT(j,3))*IV(j,3) +  b*IU(j-no,3) + (c + b*w*tau - b*w*TT(j,3))*IV(j-no,3) - (b + c)*exp(lamda*TT(j,3) - tau*lamda);
     IV(j,4)=  UTG1 - w*VTG1*TT(j,1) + (b14*D1 + b24*D2 + b34*D3 + b44*D4)*h(i) - b*IU(j-no,4) - (c + b*w*tau - b*w*TT(j,4))*IV(j-no,4) + (b + c)*exp(lamda*TT(j,4) - tau*lamda);
     IU(j,4)= (1 + w*TT(j,4))*IV(j,4) +  b*IU(j-no,4) + (c + b*w*tau - b*w*TT(j,4))*IV(j-no,4) - (b + c)*exp(lamda*TT(j,4) - tau*lamda);
     %DU(j-1,1) =2*(UTG2 - UTG1)/(h(i)); DU(j-1,2) =2*(UTG3 - UTG1)/(h(i)); DU(j-1,3) =(UTG4 - UTG1)/(h(i)); DU(j-1,4) =(U(j) - UTG1)/(b4*h(i)) - (b1*DU(j-1,1) + b2*DU(j-1,2) + b3*DU(j-1,3))/b4;
     %DV(j-1,1) =2*(VTG2 - VTG1)/(h(i)); DV(j-1,2) =2*(VTG3 - VTG1)/(h(i)); DV(j-1,3) =(VTG4 - VTG1)/(h(i)); DV(j-1,4) =(V(j) - VTG1)/(b4*h(i)) - (b1*DV(j-1,1) + b2*DV(j-1,2) + b3*DV(j-1,3))/b4;
     %== gia tri tai diem nac KT bac lien tuc==================
     IVKT(j)=  UTG1 - w*VTG1*TT(j,1) + (bkt1*D1 + bkt2*D2 + bkt3*D3 + bkt4*D4)*h(i) - b*IUKT(j-no) - (c + b*w*tau - b*w*TTKT(j))*IVKT(j-no) + (b + c)*exp(lamda*TTKT(j) - tau*lamda);
     IUKT(j)= (1 + w*TTKT(j))*IVKT(j) +  b*IUKT(j-no) + (c + b*w*tau - b*w*TTKT(j))*IVKT(j-no) - (b + c)*exp(lamda*TTKT(j) - tau*lamda);
     ErUC(j) = UCXKT(j) - IUKT(j); ErVC(j)= VCXKT(j) - IVKT(j); 
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
disp('Errors of U:'); disp(EUktc');
disp('Uniform order of U:'); disp(orderCU');
disp('Errors of V:'); disp(EVktc');
disp('Uniform order of V:'); disp(orderCV');
%Tinh sai so bang Max tren toan mien, 
disp('actual errors of U:'); disp(cu');
disp('Discrete orders of U:'); disp(orderU');
%Tinh sai so t?i diem cuoi
disp('actual errors of V:'); 
%disp(cs');
%Tinh sai so bang Max tren toan mien,  
disp(cv');
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
 
% Lap hàm dau ra: Continuous output order 3 ================
%function kq = ConO3(tm,delay,stepSZ,DH,x)
%nu = floor((tm - delay)/stepSZ) + fix(delay/stepSZ) + 3; theta = (tm - delay)/stepSZ + 3 + fix(delay/stepSZ) - nu; 
%kq = x(nu) + stepSZ*((2*(theta^3)/3 - 3*(theta^2)/2 + theta)*DH(nu,1) + ((theta^2) - 2*(theta^3)/3)*DH(nu,2) + ((theta^2) - 2*(theta^3)/3)*DH(nu,3) + (2*(theta^3)/3 - (theta^2)/2 )*DH(nu,4));
