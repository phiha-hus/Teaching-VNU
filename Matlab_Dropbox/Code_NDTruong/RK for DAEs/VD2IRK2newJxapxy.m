function VD2IRK2newJxapxy
%% tieu chuan dung: sai so giua hai phép lap liên tiep nho hon h^6 hoac so phép lap lon hon 50
clear all
format short e;
% Truong hop matran Jacobi rut gon(lay xap xi)
%h=input('Nhap buoc thoi gian h = ');
%N=8/h;
h0=input('nhap buoc thoi gian dau h0 = ');
%z =input('Nhap so lan lap Newton z = ');
m=input('Nhap so lan giai m =');
time = 1;
for k =1:m
    h(k)=2*h0/(2^k);
    N = time/h(k);
    %N = time/h;
    for i=1:N+1
        t(i)=(i-1)*h(k);
    end;
    %U xap xi cua x1, V xap xi cua x2 
    U(1)=1;
    V(1)=0; UCX(1)=1; VCX(1)=0;
    for i=1:N
         ttg(1)=t(i)+ h(k)/3; ttg(2)=t(i)+ h(k);
         BTG(1,1)= U(i);  BTG(1,2)= t(i)*U(i); BTG(2,1)= exp(-t(i)); BTG(2,2)= -1;
         TGB = inv(BTG); 
         A(1,1)=5*TGB(1,1)/12; A(1,2)=5*TGB(1,2)/12; A(1,3)= -TGB(1,1)/12; A(1,4)= -TGB(1,2)/12;
         A(2,1)=5*TGB(2,1)/12; A(2,2)=5*TGB(2,2)/12; A(2,3)= -TGB(2,1)/12; A(2,4)= -TGB(2,2)/12;
         A(3,1)=3*TGB(1,1)/4; A(3,2)=3*TGB(1,2)/4; A(3,3)=TGB(1,1)/4; A(3,4)=TGB(1,2)/4;
         A(4,1)=3*TGB(2,1)/4; A(4,2)=3*TGB(2,2)/4; A(4,3)=TGB(2,1)/4; A(4,4)=TGB(2,2)/4;
         UVTG(1)=U(i); UVTG(2)=V(i);UVTG(3)=U(i); UVTG(4)=V(i); K1=0; K2=0;
         %for e=1:z
         dem = 0; maxci=1;
         while dem <=50
           if maxci >= h(k)^6
           dem = dem +1;
         F(1)= fu1(h(k),ttg(1),UVTG(1),UVTG(2),U(i),V(i),K1);
         F(2)= fu2(ttg(1),UVTG(1),UVTG(2));
         F(3)= fu1(h(k),ttg(2),UVTG(3),UVTG(4),U(i),V(i),K2);
         F(4)= fu2(ttg(2),UVTG(3),UVTG(4));
         %F(2)=(3*sin(ttg(1))+ sin(ttg(2)))/2 - 2*sin(t(i)); F(4)=(-9*sin(ttg(1))+ 5*sin(ttg(2)))/2 + 2*sin(t(i));        
         % phép lap Newton 
         UVC= UVTG' - A*(F');
         for j=1:4
         Hieu(j) = abs(UVC(j)-UVTG(j));
         end
         maxci = max(Hieu);
         UVTG = UVC'; 
         K1 = (3*(UVTG(1)+ ttg(1)*UVTG(2) - U(i)- t(i)*V(i))+ (UVTG(3) + ttg(2)*UVTG(4) - U(i)- t(i)*V(i)))/(2*h(k));
         K2 =(-9*(UVTG(1)+ ttg(1)*UVTG(2) - U(i)- t(i)*V(i))+ 5*(UVTG(3)+ ttg(2)*UVTG(4) - U(i)- t(i)*V(i)))/(2*h(k));
            else
            break
           end
         end
         %UVD(1)= U(i); UVD(2)= V(i); UVD(3)= U(i); UVD(4)= V(i); 
        
         U(i+1)= (U(i)+ V(i)*t(i) + (3*K1 + K2)*h(k)/4 - t(i+1)*sin(t(i+1))+ t(i+1))/(1+ t(i+1)*exp(-t(i+1)));
         V(i+1)= U(i+1)*exp(-t(i+1)) + sin(t(i+1))- 1;
         UCX(i+1)=exp(t(i+1));    VCX(i+1)=sin(t(i+1));
         SSU(i+1)=abs(U(i+1)-UCX(i+1));    SSV(i+1)=abs(V(i+1)-VCX(i+1));
         %UTG(0)=U(i); VTG(0)=V(i);
         % A là ma tr?n Jacobian
    end;
    a(k) = max(SSU);
    b(k) = max(SSV);
end
for e=1:(m-1)
   orderU(e)=log(a(e)/a(e+1))/log(2);
   orderV(e)=log(b(e)/b(e+1))/log(2);
end;
disp('actual errors of U:'); disp(a');
disp('error orders of U:'); disp(orderU');
disp('actual errors of V:'); disp(b');
disp('error orders of V:'); disp(orderV');
function D=fu1(h,t,U,V,u0,v0,K)
D=U*K*h - U*V*h - h*U*V*exp(t) - h*(exp(2*t)+ t*cos(t)*exp(t) - sin(t)*exp(2*t));
function E=fu2(t,U,V)
E = U*exp(-t) - V + sin(t)- 1;