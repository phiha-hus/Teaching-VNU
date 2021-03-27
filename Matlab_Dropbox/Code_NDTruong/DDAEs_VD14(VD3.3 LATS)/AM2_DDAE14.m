function AM2_DDAE14
clear all
format short e;
beta0 = 5/12; beta1 = 2/3; beta2 = -1/12;
h0= input('nhap buoc thoi gian dau (la uoc cua Pi), h0 = ');
m=7; %input('Nhap so lan giai m =');
time = 10*pi;
tau = pi;
for i =1:m
    h(i)=2*h0/(2^i);
 % tong so diem chia la
 M = fix(time/h(i)); 
 % So diem chia trên mot khoang tre ====================
 no = fix(tau/h(i));
 % chia luoi thoi gia va gan gia tri chinh xac ============================
    for k = 1:M+no+3
        t(k)=(k-no -3)*h(i);
        UCX(k)= exp(-t(k));   VCX(k)= sin(t(k));
    end;
 % gan gia tri ban dau ===================================================   
    for k=1:no+3
        U(k)= UCX(k);    V(k)=VCX(k); 
        SSU(k)= 0; SSV(k)= 0; 
    end;
 %========= GT ban dau cua Dao Ham EX====================
    D(no+2)= -exp(-t(no+2)) + 2*t(no +2)*sin(t(no+2)) + t(no+2)*t(no +2)*cos(t(no+2)) + 2*sin(2*t(no+2)); 
    D(no+3)= -exp(-t(no+3)) + 2*t(no +3)*sin(t(no+3)) + t(no+3)*t(no +3)*cos(t(no+3)) + 2*sin(2*t(no+3)); 
    for j=no+3:M+no+2
        TG(1)= U(j); TG(2)= V(j);
        % th? tuc l?p Neuton
        dem = 0; maxci=1;
         while dem <=10
           if maxci >= (h(i))^4
           dem = dem +1;
           TGD = TG;
         %UVD(1)= U(i); UVD(2)= V(i); UVD(3)= U(i); UVD(4)= V(i); 
         A(1,1) = 2*TGD(1) + TGD(2)*(t(j+1)*t(j+1) + 2*sin(t(j+1)) - h(i)*beta0*(2*t(j+1) + 2*cos(t(j+1))+ exp(-t(j+1)))) - U(j) - (t(j)*t(j) + 2*sin(t(j)))*V(j) - (beta0*sin(2*t(j+1)) + beta1*D(j) + beta2*D(j-1))*h(i); 
         A(1,2) = TGD(1)*(t(j+1)*t(j+1) + 2*sin(t(j+1)) - h(i)*beta0*(2*t(j+1) + 2*cos(t(j+1)) + exp(-t(j+1))));
         A(2,1) = exp(t(j+1)); A(2,2)= -1; 
         %JND = inv(A);
         F(1)= FA1(TGD(1),U(j),TGD(2),V(j),noisuy1(4,tau,t(j+1),h(i),t,V),D(j),D(j-1),beta0,beta1,beta2,t(j+1),h(i));
         F(2)= FA2(TGD(1),TGD(2),noisuy1(4,tau,t(j+1),h(i),t,V),t(j+1));
         %TG0 = TGD' - JND*(F');
         TH=A\(F');
         TG = TGD - TH';
         %TG = TG0;
         for k=1:2
         Hieu(k) = abs(TG(k)-TGD(k));
         end
         maxci = max(Hieu);
            else
                break
            end
         end
        U(j+1) = TG(1); V(j+1) = TG(2);
        D(j+1)= (U(j+1) + (t(j+1)*t(j+1) + 2*sin(t(j+1)))*V(j+1) - U(j) - (t(j)*t(j) + 2*sin(t(j)))*V(j))/(beta0*h(i)) - (beta1*D(j) + beta2*D(j-1))/beta0;      
        SSU(j+1)=abs(U(j+1) - UCX(j+1));    SSV(j+1)= abs(V(j+1) - VCX(j+1));
    end;
    a(i)=max(SSU);
    b(i)=max(SSV);
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
function D=FA1(x1,x2,y1,y2,qky,dh1,dh2,bt0,bt1,bt2,t,h1)
D=x1*x1 + x1*y1*(t*t + 2*sin(t) - h1*bt0*(2*t + 2*cos(t) + exp(-t))) - x1*(x2 + y2*((t-h1)*(t - h1) + 2*sin(t-h1)) + h1*(bt0*sin(2*t) + bt1*dh1 + bt2*dh2)) - h1*bt0*qky*exp(-2*t) - h1*bt0*exp(-t)*(t*t*cos(t) - exp(-t));
function E=FA2(x1,y1,qky,t)
E = x1*exp(t) - y1 - qky - 1;
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