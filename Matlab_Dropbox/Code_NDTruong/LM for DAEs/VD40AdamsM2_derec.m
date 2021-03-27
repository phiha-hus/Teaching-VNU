function VD40AdamsM2_derec
clear all
format short e;
%beta0 = 5/12; beta1 = 2/3; beta2 = -1/12; % Adams Moulton 2 Method
beta0 = 1/2; beta1 = 1/2;  % Trapezoidal Method
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
    U(1)=1; V(1)=0; %U(2)= UCX(2); V(2)= VCX(2);% U(3)= UCX(3); V(3)= VCX(3); 
    SSU(1)=0; SSV(1)=0; %SSU(2)= 0; SSV(2)= 0; %SSU(3)= 0; SSV(3)= 0;
    DU(1)= -1; DV(1)= 1;% DU(2)= - exp(-t(2)); DV(2)= cos(t(2)); 
    for i=1:N
        %UTG1=U(i); VTG1=V(i);
        TG(1)= U(i); TG(2)= V(i);
        % th? tuc l?p Neuton
        dem = 0; maxci=1;
         while dem <=10
           if maxci >= (h(k))^3
           dem = dem +1;
           TGD = TG;
         %----- Trapezoidal method
          A(1,1) = 2*TGD(1) + TGD(2)*(t(i+1)*t(i+1) + 2*sin(t(i+1))- h(k)*beta0*(exp(-t(i+1)))) - U(i) - (t(i+1)*t(i+1) + 2*sin(t(i+1)))*V(i) - (beta0*sin(2*t(i+1)) + beta1*DU(i))*h(k) - (t(i+1)*t(i+1) + 2*sin(t(i+1)))*(beta1*DV(i))*h(k); 
         %------- Adams Moulton 2 method
         %A(1,1) = 2*TGD(1) + TGD(2)*(t(i+1)*t(i+1) + 2*sin(t(i+1))- h(k)*beta0*(exp(-t(i+1)))) - U(i) - (t(i+1)*t(i+1) + 2*sin(t(i+1)))*V(i) - (beta0*sin(2*t(i+1)) + beta1*DU(i) + beta2*DU(i-1))*h(k) - (t(i+1)*t(i+1) + 2*sin(t(i+1)))*(beta1*DV(i) + beta2*DV(i-1))*h(k); 
         A(1,2) = TGD(1)*(t(i+1)*t(i+1) + 2*sin(t(i+1)) - h(k)*beta0*(exp(-t(i+1)))) + h(k)*beta0*(exp(-2*t(i+1)));
         A(2,1) = exp(t(i+1)); A(2,2)= -1; 
         %JND = inv(A);
         %----- Trapezoidal method
         F(1)= FA1Tra(TGD(1),U(i),TGD(2),V(i),DU(i),DV(i),beta0,beta1,t(i+1),h(k));
         %------- Adams Moulton 2 method
         %F(1)= FA1Adm2(TGD(1),U(i),TGD(2),V(i),DU(i),DU(i-1),DV(i),DV(i-1),beta0,beta1,beta2,t(i+1),h(k));
         F(2)= FA2(TGD(1),TGD(2),t(i+1));
         %TG0 = TGD' - JND*(F');
         TH=A\(F');
         TG = TGD - TH';
         %TG = TG0;
         for j=1:2
         Hieu(j) = abs(TG(j) - TGD(j));
         end
         maxci = max(Hieu);
            else
                break
            end
         end
        U(i+1) = TG(1); V(i+1) = TG(2);
        %------- Trapezoidal method
        DU(i+1)= (U(i+1) - U(i))/(beta0*h(k)) - (beta1*DU(i))/beta0; 
        DV(i+1)= (V(i+1) - V(i))/(beta0*h(k)) - (beta1*DV(i))/beta0;
        %------- Adams Moulton 2 method
        %DU(i+1)= (U(i+1) - U(i))/(beta0*h(k)) - (beta1*DU(i) + beta2*DU(i-1))/beta0; 
        %DV(i+1)= (V(i+1) - V(i))/(beta0*h(k)) - (beta1*DV(i) + beta2*DV(i-1))/beta0;
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
%--------- Trapezoidal Method
function D=FA1Tra(x1,x2,y1,y2,du1,dv1,bt0,bt1,t,h1)
D=x1*x1 + x1*y1*(t*t + 2*sin(t) - h1*bt0*(exp(-t))) - x1*(x2 + y2*(t*t + 2*sin(t)) + h1*(bt1*du1) + h1*(t*t + 2*sin(t))*(bt1*dv1)) - h1*bt0*(t*t*cos(t) - exp(-t))*exp(-t) + h1*bt0*(y1*exp(-2*t) - x1*sin(2*t));
%-------- Adams Moulton 2-step method
function D=FA1Adm2(x1,x2,y1,y2,du1,du2,dv1,dv2,bt0,bt1,bt2,t,h1)
D=x1*x1 + x1*y1*(t*t + 2*sin(t) - h1*bt0*(exp(-t))) - x1*(x2 + y2*(t*t + 2*sin(t)) + h1*(bt1*du1 + bt2*du2) + h1*(t*t + 2*sin(t))*(bt1*dv1 + bt2*dv2)) - h1*bt0*(t*t*cos(t) - exp(-t))*exp(-t) + h1*bt0*(y1*exp(-2*t) - x1*sin(2*t));
function E=FA2(x1,y1,t)
E = x1*exp(t) - y1 + sin(t) - 1;