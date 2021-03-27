function VD14IEULERWOR01
clear all
format short e;
%h=input('Nhap buoc thoi gian h = ');
%N=8/h;
h0=input('nhap buoc thoi gian dau h0 = ');
z =input('Nhap so lan lap Newton z = ');
m =input('Nhap so lan giai m = ');
time = 2;
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
         %ttg(1)=t(i)+ h/3; ttg(2)=t(i)+ h;
         U(i+1)=U(i); V(i+1)=V(i);
         for j = 1:z
         UTG=U(i+1); VTG = V(i+1);
         % A là ma tr?n Jacobian
         A(1,1)=2*UTG + t(i+1)*VTG - U(i) - t(i)*V(i)- h(k)*VTG*exp(t(i+1)) - h(k)*VTG; A(1,2)= UTG*t(i+1)- h(k)*UTG - h(k)*UTG*exp(t(i+1));
         A(2,1) = exp(-t(i+1)); A(2,2)= -1; 
         % JND là ma tr?n ngh?ch ð?o Jacobian
         JND = inv(A); 
         d(1)= fun01(h(k),t(i+1),UTG,VTG,U(i),V(i)); 
         d(2)= fun02(t(i+1),UTG,VTG);
         UVD(1)= UTG; UVD(2)= VTG; 
         % phép lap Newton 
         UVC= UVD' - JND*(d');
         %UVC2= UVC1 - A*(F');
         %UVC= UVC';
         U(i+1)= UVC(1); V(i+1)= UVC(2);
         end;
         UCX(i+1)= exp(t(i+1));    VCX(i+1)=sin(t(i+1));
         SSU(i+1)= abs(U(i+1)- UCX(i+1));    SSV(i+1)=abs(V(i+1)- VCX(i+1));
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
function D=fun01(h,t,U,V,u0,v0)
D=U*U + t*U*V -U*u0 - U*v0*(t-h)- h*U*V - h*U*V*exp(t) - h*exp(2*t)- h*t*cos(t)*exp(t) + h*sin(t)*exp(2*t);
function B=fun02(t,U,V)
B = U*exp(-t) - V + sin(t)-1;