function VD2IRK2new
%% tieu chuan dung: sai so giua hai ph�p lap li�n tiep nho hon h^6 hoac so ph�p lap lon hon 10
clear all
format short e;
%h=input('Nhap buoc thoi gian h = ');
%N=8/h;
h0=input('nhap buoc thoi gian dau h0 = ');
m=input('Nhap so lan giai m =');
%z =input('Nhap so lan lap Newton z = ');
%===So buoc l?p Newton 6,...10.....
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
         UVTG(1)=U(i); UVTG(2)=V(i);UVTG(3)=U(i); UVTG(4)=V(i);K1=0; K2=0;
         % ph�p lap Newton 
         %for j = 1:z
         dem = 0; maxci=1;
         while dem <=10
           if maxci >= h(k)^6
           dem = dem +1;    
         %UVD(1)= U(i); UVD(2)= V(i); UVD(3)= U(i); UVD(4)= V(i); 
         A(1,1)=3*UVTG(1)+ 3*UVTG(2)*ttg(1)/2 + UVTG(3)/2 + UVTG(4)*ttg(2)-2*U(i)- 2*V(i)*t(i)- h(k)*UVTG(2)*exp(ttg(1)) - h(k)*UVTG(2); 
         A(1,2)=3*UVTG(1)*ttg(1)/2 - h(k)*UVTG(1)*exp(ttg(1))- h(k)*UVTG(1); A(1,3)= UVTG(1)/2; A(1,4)= UVTG(1)*ttg(2)/2;
         A(2,1) = exp(-ttg(1)); A(2,2)= -1; A(2,3)= 0; A(2,4)= 0;
         A(3,1)= -9*UVTG(3)/2; A(3,2)= -9*UVTG(3)*ttg(1)/2; 
         A(3,3)= -9*UVTG(1)/2 - 9*UVTG(2)*ttg(1)/2 + 5*UVTG(3)+ 5*UVTG(4)*ttg(2)/2 + 2*U(i)+ 2*V(i)*t(i)- h(k)*UVTG(4)*exp(ttg(2))- h(k)*UVTG(4); 
         A(3,4)=5*UVTG(3)*ttg(2)/2 - h(k)*UVTG(3)*exp(ttg(2))- h(k)*UVTG(3);
         A(4,1)=0; A(4,2)= 0; A(4,3)= exp(-ttg(2));A(4,4)= -1;
         JND = inv(A);
         F(1)= UVTG(1)*h(k)*K1 - UVTG(1)*UVTG(2)*h(k) - UVTG(1)*UVTG(2)*h(k)*exp(ttg(1)) - fu3(h(k),ttg(1));
         F(2)= fu2(ttg(1),UVTG(1),UVTG(2));
         F(3)= UVTG(3)*h(k)*K2 - UVTG(3)*UVTG(4)*h(k) - UVTG(3)*UVTG(4)*h(k)*exp(ttg(2)) - fu3(h(k),ttg(2));
         F(4)= fu2(ttg(2),UVTG(3),UVTG(4));
         %F(2)=(3*sin(ttg(1))+ sin(ttg(2)))/2 - 2*sin(t(i)); F(4)=(-9*sin(ttg(1))+ 5*sin(ttg(2)))/2 + 2*sin(t(i));        
         UVC= UVTG' - JND*(F');
         UVTG = UVC';
         K1 = (3*(UVTG(1)+ ttg(1)*UVTG(2) - U(i)- t(i)*V(i))+ (UVTG(3) + ttg(2)*UVTG(4) - U(i)- t(i)*V(i)))/(2*h(k));
         K2 =(-9*(UVTG(1)+ ttg(1)*UVTG(2) - U(i)- t(i)*V(i))+ 5*(UVTG(3)+ ttg(2)*UVTG(4) - U(i)- t(i)*V(i)))/(2*h(k));
            else
                break
            end
         end
         %end
         U(i+1)= (U(i)+ V(i)*t(i) + (3*K1 + K2)*h(k)/4 - t(i+1)*sin(t(i+1))+ t(i+1))/(1+ t(i+1)*exp(-t(i+1)));
         V(i+1)= U(i+1)*exp(-t(i+1)) + sin(t(i+1))- 1;
         UCX(i+1)=exp(t(i+1));    VCX(i+1)=sin(t(i+1));
         SSU(i+1)=abs(U(i+1)-UCX(i+1));    SSV(i+1)=abs(V(i+1)-VCX(i+1));        
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
%subplot(2,2,3);
%plot(t(1:(N+1)),U(1:(N+1)),'+',t(1:(N+1)),UCX(1:(N+1)),'-')
%title ('nghiem xap xi(+),nghiem dung (-), Ham e')
%subplot(2,2,4);
%plot(t(1:(N+1)),V(1:(N+1)),'+',t(1:(N+1)),VCX(1:(N+1)),'-')
%title ('nghiem xap xi(+),nghiem dung (-), Ham sin')
function D=fu3(h,t)
D=h*(exp(2*t)+ t*cos(t)*exp(t) - sin(t)*exp(2*t));
function E=fu2(t,U,V)
E = U*exp(-t) - V + sin(t)- 1;
