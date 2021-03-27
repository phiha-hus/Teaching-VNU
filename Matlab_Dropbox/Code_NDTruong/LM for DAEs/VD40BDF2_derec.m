function VD40BDF2_derec
clear all
format short e;
beta0 = 2/3; alpha1 = - 4/3; alpha2 = 1/3;%input('Nhap alpha = ');
h0=0.1; %input('nhap buoc thoi gian dau h0 = ');
m=7; %input('Nhap so lan giai m =');
time = 100;
for k =1:m
    h(k)=2*h0/(2^k);
    N = time/h(k);
    for i=1:N+1
        t(i)=(i-1)*h(k);  UCX(i)=exp(-t(i)); VCX(i)=sin(t(i));
    end;
    %U xap xi cua x1, V xap xi cua x2 
    U(1)=1; V(1)=0; U(2)= UCX(2); V(2)= VCX(2);% U(3)= UCX(3); V(3)= VCX(3); 
    SSU(1)=0; SSV(1)=0; SSU(2)= 0; SSV(2)= 0; %SSU(3)= 0; SSV(3)= 0;
    %DU(1)= -1; DV(1)= 1;% DU(2)= - exp(-t(2)); DV(2)= cos(t(2)); 
    for i=2:N
        %UTG1=U(i); VTG1=V(i);
        TG(1)= U(i); TG(2)= V(i);
        % th? tuc l?p Neuton
        dem = 0; maxci=1;
         while dem <=10
           if maxci >= (h(k))^3
           dem = dem +1;
           TGD = TG;
         %UVD(1)= U(i); UVD(2)= V(i); UVD(3)= U(i); UVD(4)= V(i); 
         A(1,1) = 2*TGD(1) + TGD(2)*(t(i+1)*t(i+1) + 2*sin(t(i+1))- h(k)*beta0*(exp(-t(i+1)))) + alpha1*U(i) + alpha2*U(i-1) + (t(i+1)*t(i+1) + 2*sin(t(i+1)))*(alpha1*V(i) + alpha2*V(i-1)) - h(k)*beta0*(sin(2*t(i+1))); 
         A(1,2) = TGD(1)*(t(i+1)*t(i+1) + 2*sin(t(i+1)) - h(k)*beta0*(exp(-t(i+1)))) + h(k)*beta0*(exp(-2*t(i+1)));
         A(2,1) = exp(t(i+1)); A(2,2)= -1; 
         %JND = inv(A);
         F(1)= FA1(TGD(1),U(i),U(i-1),TGD(2),V(i),V(i-1),beta0,alpha1,alpha2,t(i+1),h(k));
         F(2)= FA2(TGD(1),TGD(2),t(i+1));
         %TG0 = TGD' - JND*(F');
         TH=A\(F');
         TG = TGD - TH';
         %TG = TG0;
         for j=1:2
         Hieu(j) = abs(TG(j)-TGD(j));
         end
         maxci = max(Hieu);
            else
                break
            end
         end
        U(i+1) = TG(1); V(i+1) = TG(2);
        %DU(i+1)= (U(i+1) - U(i))/(beta0*h(k)) - (beta1*DU(i))/beta0; 
        %DV(i+1)= (V(i+1) - V(i))/(beta0*h(k)) - (beta1*DV(i))/beta0;
        SSU(i+1)=abs(U(i+1) - UCX(i+1));    SSV(i+1)= abs(V(i+1) - VCX(i+1));
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
disp('actual errors of U:'); disp(a');
disp('error orders of U:'); disp(orderCXU');
disp('actual errors of V:'); disp(b');
disp('error orders of V:'); disp(orderCXV');
function D=FA1(x1,x2,x3,y1,y2,y3,bt0,alp1,alp2,t,h1)
D=x1*x1 + x1*y1*(t*t + 2*sin(t) - h1*bt0*(exp(-t))) + x1*(alp1*x2 + alp2*x3 + (alp1*y2 + alp2*y3)*(t*t + 2*sin(t)) - h1*bt0*sin(2*t)) + h1*bt0*y1*exp(-2*t)  - h1*bt0*(t*t*cos(t) - exp(-t))*exp(-t);
function E=FA2(x1,y1,t)
E = x1*exp(t) - y1 + sin(t) - 1;