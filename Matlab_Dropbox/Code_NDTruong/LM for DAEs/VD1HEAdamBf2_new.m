function VD1HEAdamBf2_new
clear all
format short e;
lamda =input('Nhap lamda = ');
%N=8/h;
w = input('Nhap w = ');
beta1 = 1.5; beta2 = -0.5;%input('Nhap alpha = ');
h0=0.1; %input('nhap buoc thoi gian dau h0 = ');
m=6; %input('Nhap so lan giai m =');
time = 5;
for k =1:m
    h(k)=2*h0/(2^k);
    N = time/h(k);
    for i=1:N+1
        t(i)=(i-1)*h(k);  VCX(i)=exp(lamda*t(i)); UCX(i)=(1 + w*t(i))*exp(lamda*t(i));
    end;
    %U xap xi cua x1, V xap xi cua x2 
    U(1)=1;  V(1)=1; SSU(1)=0; SSV(1)=0; DH(1)=lamda*exp(-lamda*h(k));
    for i=1:N
        %UTG1=U(i); VTG1=V(i);
        %VP(i) = 
        V(i+1) = (1 + beta1*lamda*h(k))*U(i) - (w*t(i))*(1 + beta1*lamda*h(k))*V(i) + h(k)*beta2*DH(i);
        U(i+1) = (1 + w*t(i+1))*V(i+1); 
        DH(i+1)= (U(i+1) - w*t(i+1)*V(i+1) - U(i) + w*t(i)*V(i))/(beta1*h(k))- (beta2*DH(i))/beta1;      
        SSU(i+1)=abs(U(i+1)-UCX(i+1));    SSV(i+1)=abs(V(i+1)-VCX(i+1));
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