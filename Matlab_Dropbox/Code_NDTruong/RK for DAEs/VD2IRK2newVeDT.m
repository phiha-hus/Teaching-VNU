function VD2IRK2newVeDT
%%%% Ve Do thi.................
clear all
format short e;
h=input('Nhap buoc thoi gian h = ');
z =input('Nhap so lan lap Newton z = ');
%h0=input('nhap buoc thoi gian dau h0 = ');
%m=input('Nhap so lan giai m =');
time = 2;
%for k =1:m
    %h(k)=2*h0/(2^k);
    N = time/h;
    for i=1:N+1
        t(i)=(i-1)*h;
    end;
    %U xap xi cua x1, V xap xi cua x2 
    U(1)=1;
    V(1)=0;UCX(1)=1;VCX(1)=0;
    for i=1:N
         ttg(1)=t(i)+ h/3; ttg(2)=t(i)+ h;
         %UTG(0)=U(i); VTG(0)=V(i);
         BTG(1,1)= U(i);  BTG(1,2)= t(i)*U(i); BTG(2,1)= exp(-t(i)); BTG(2,2)= -1;
         TGB = inv(BTG); 
         A(1,1)=5*TGB(1,1)/12; A(1,2)=5*TGB(1,2)/12; A(1,3)= -TGB(1,1)/12; A(1,4)= -TGB(1,2)/12;
         A(2,1)=5*TGB(2,1)/12; A(2,2)=5*TGB(2,2)/12; A(2,3)= -TGB(2,1)/12; A(2,4)= -TGB(2,2)/12;
         A(3,1)=3*TGB(1,1)/4; A(3,2)=3*TGB(1,2)/4; A(3,3)=TGB(1,1)/4; A(3,4)=TGB(1,2)/4;
         A(4,1)=3*TGB(2,1)/4; A(4,2)=3*TGB(2,2)/4; A(4,3)=TGB(2,1)/4; A(4,4)=TGB(2,2)/4;
         UVTG(1)=U(i); UVTG(2)=V(i);UVTG(3)=U(i); UVTG(4)=V(i);K1=0; K2=0;
         for j = 1:z
         %UVD(1)= U(i); UVD(2)= V(i); UVD(3)= U(i); UVD(4)= V(i); 
         F(1)= fu1(h,ttg(1),UVTG(1),UVTG(2),U(i),V(i),K1);
         F(2)= fu2(ttg(1),UVTG(1),UVTG(2));
         F(3)= fu1(h,ttg(2),UVTG(3),UVTG(4),U(i),V(i),K2);
         F(4)= fu2(ttg(2),UVTG(3),UVTG(4));
         %F(2)=(3*sin(ttg(1))+ sin(ttg(2)))/2 - 2*sin(t(i)); F(4)=(-9*sin(ttg(1))+ 5*sin(ttg(2)))/2 + 2*sin(t(i));        
         % phép lap Newton 1 lan
         UVC= UVTG' - A*(F');
         UVTG = UVC';
         K1 = (3*(UVTG(1)+ ttg(1)*UVTG(2) - U(i)- t(i)*V(i))+ (UVTG(3) + ttg(2)*UVTG(4) - U(i)- t(i)*V(i)))/(2*h);
         K2 =(-9*(UVTG(1)+ ttg(1)*UVTG(2) - U(i)- t(i)*V(i))+ 5*(UVTG(3)+ ttg(2)*UVTG(4) - U(i)- t(i)*V(i)))/(2*h);
         end
         U(i+1)= (U(i)+ V(i)*t(i) + (3*K1 + K2)*h/4 - t(i+1)*sin(t(i+1))+ t(i+1))/(1+ t(i+1)*exp(-t(i+1)));
         V(i+1)= U(i+1)*exp(-t(i+1)) + sin(t(i+1))- 1;
         UCX(i+1)=exp(t(i+1));    VCX(i+1)=sin(t(i+1));
        SSU(i+1)=abs(U(i+1)-UCX(i+1));    SSV(i+1)=abs(V(i+1)-VCX(i+1));
    end;
    a = max(SSU)
    b = max(SSV)
%end
%for e=1:(m-1)
    %orderCXU(e)=log(a(e)/a(e+1))/log(2);
   % orderCXV(e)=log(b(e)/b(e+1))/log(2);
% end;
%c=a';
%d=order';
%D=[b; c]
%disp('actual errors of U:'); disp(a');
%disp('error orders of U:'); disp(orderCXU');
%disp('actual errors of V:'); disp(b');
%disp('error orders of V:'); disp(orderCXV');
subplot(2,2,3);
plot(t(1:(N+1)),U(1:(N+1)),'+',t(1:(N+1)),UCX(1:(N+1)),'-')
title ('nghiem xap xi(+),nghiem dung (-), Ham e')
subplot(2,2,4);
plot(t(1:(N+1)),V(1:(N+1)),'+',t(1:(N+1)),VCX(1:(N+1)),'-')
title ('nghiem xap xi(+),nghiem dung (-), Ham sin')
function D=fu1(h,t,U,V,u0,v0,K)
D=U*K*h - U*V*h - h*U*V*exp(t) - h*(exp(2*t)+ t*cos(t)*exp(t) - sin(t)*exp(2*t));
function E=fu2(t,U,V)
E = U*exp(-t) - V + sin(t)- 1;