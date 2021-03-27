function VD1HELM3_derec
clear all
format short e;
lamda =input('Nhap lamda = ');
%N=8/h;
w = input('Nhap w = ');
beta1 = 0.5; beta2 = 1.5; beta3 = -1;%input('Nhap alpha = ');
h0=0.1; %input('nhap buoc thoi gian dau h0 = ');
m=8; %input('Nhap so lan giai m =');
time = 5;
for k =1:m
    h(k)=2*h0/(2^k);
    N = time/h(k);
    for i=1:N+1
        t(i)=(i-1)*h(k);  VCX(i)=exp(lamda*t(i)); UCX(i)=(1 + w*t(i))*exp(lamda*t(i));
    end;
    % Giá tr? kh?i t?o cho U, V, DH, SS.... 
    U(1)=1; U(2)= UCX(2); U(3)= UCX(3);  V(1)=1; V(2)= VCX(2); V(3)= VCX(3); 
    SSU(1)=0; SSU(2)= 0; SSU(3)= 0; SSV(1)=0; SSV(2)= 0; SSV(3)= 0;
    DHV(1)=lamda; DHV(2)= lamda*exp(lamda*t(2));
    DHU(1)=lamda + w; DHU(2)= (lamda + w + w*t(2))*exp(lamda*t(2));
    for i=3:N
        VP(i) = U(i) - w*t(i)*V(i) + h(k)*beta2*DHU(i-1) + h(k)*beta3*DHU(i-2) - h(k)*beta2*t(i)*w*DHV(i-1)- h(k)*beta3*t(i)*w*DHV(i-2)+ beta1*h(k)*lamda*U(i)+ beta1*h(k)*(w - w*lamda*t(i))*V(i);
        V(i+1) = VP(i)/(1 + w*h(k));% U(i) - w*t(i)*V(i) + beta2*h(k)*DH(i-1) + beta3*h(k)*DH(i-2) + h(k)*beta1*lamda*U(i) - h(k)*beta1*lamda*w*t(i)*V(i);
        U(i+1) = (1 + w*t(i+1))*V(i+1); 
        DHU(i)= (U(i+1) - U(i))/(beta1*h(k)) - (beta2*DHU(i-1))/beta1 - (beta3*DHU(i-2))/beta1;
        DHV(i)= (V(i+1) - V(i))/(beta1*h(k)) - (beta2*DHV(i-1))/beta1 - (beta3*DHV(i-2))/beta1;
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