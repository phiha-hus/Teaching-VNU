function VD1HERK2stand
clear all
format short e;
lamda = -1;% input('Nhap lamda = ');
%N=8/h;
w = input('Nhap w = ');
alpha = 1; %input('Nhap alpha = ');
h0=0.1; %input('nhap buoc thoi gian dau h0 = ');
m=8; %input('Nhap so lan giai m =');
time = 5;
for k =1:m
    h(k)=2*h0/(2^k);
    N = time/h(k);
    for i=1:N+1
        t(i)=(i-1)*h(k);
        VCX(i)=exp(lamda*t(i));    UCX(i)=(1 + w*t(i))*VCX(i);
    end;
    %U xap xi cua x1, V xap xi cua x2 
    U(1)=1;   V(1)=1; SSU(1)=0; SSV(1)=0;
    for i=1:N
        ttg1=t(i);ttg2=t(i)+ alpha*h(k);
        UTG1=U(i); VTG1=V(i);
        VTG2 = (UTG1 - ttg1*w*VTG1 + h(k)*alpha*(lamda*UTG1 + w*VTG1*(1 - lamda*ttg1)))/(1 + w*alpha*h(k));
        UTG2 = (1 + w*ttg2)*VTG2;
        A1 = (UTG2 - UTG1)/(h(k)*alpha);
        B1 = (VTG2 - VTG1)/(h(k)*alpha);
        V(i+1) = (U(i) - w*ttg2*V(i) + h(k)*lamda*(UTG2 - w*ttg2*VTG2)/(2*alpha)+ h(k)*w*VTG2/(2*alpha) - h(k)*A1*(1 - 2*alpha)/(2*alpha) + w*ttg2*h(k)*B1*(1 - 2*alpha)/(2*alpha))/(1 + w*(t(i+1)-ttg2));
        U(i+1) = (1 + w*t(i+1))*V(i+1);        
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