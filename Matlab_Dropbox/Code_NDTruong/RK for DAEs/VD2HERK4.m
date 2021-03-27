function VD2HERK4
clear all
format short e;
%h=input('Nhap bý?c thoi gian h = ');
%N=8/h;
h0=input('nhap buoc thoi gian dau h0 = ');
m=input('Nhap so lan giai m =');
%batdau = clock;
time = 1;
for k =1:m
    h(k)=2*h0/(2^k);
    N = time/h(k);
    for i=1:N+1
        t(i)=(i-1)*h(k);
    end;
    %U xap xi cua x1, V xap xi cua x2 
    U(1)=1;
    V(1)=0; UCX(1)=1; VCX(1)=0; SSU(1)=0; SSV(1)=0;
    for i=1:N
        ttg(1)=t(i);ttg(2)=t(i)+ h(k)/2; ttg(3)=ttg(2); ttg(4)=t(i)+h(k);
        UTG(1)=U(i); VTG(1)=V(i);
        for j=2:3
            VP1=h(k)*funt2(ttg(j-1),UTG(j-1),VTG(j-1))/(2*UTG(j-1))+ UTG(1)+ ttg(j-1)*VTG(1);
            VP2=1- sin(ttg(j));
            MT(1,1)=1; MT(2,2)=-1;MT(1,2)=ttg(j-1);MT(2,1)= exp(-ttg(j));
            ABM=[UTG(j-1);VTG(j-1)]- inv(MT)*(MT*[UTG(j-1);VTG(j-1)]-[VP1;VP2]);
            UTG(j)=ABM(1,1);   VTG(j)=ABM(2,1);
            DHU(j-1)=2*(UTG(j)-UTG(1))/h(k);DHV(j-1)=2*(VTG(j)-VTG(1))/h(k);
        end;
        %buoc thu 4
        VP3= h(k)*funt2(ttg(3),UTG(3),VTG(3))/UTG(3)+ UTG(1)+ ttg(3)*VTG(1);
        VP4=1- sin(ttg(4));
        MT(1,1)=1; MT(2,2)=-1;MT(1,2)= ttg(3);MT(2,1)= exp(-ttg(4));
        AB4=[UTG(3);VTG(3)]- inv(MT)*(MT*[UTG(3);VTG(3)]-[VP3;VP4]);
        UTG(4)=AB4(1,1);VTG(4)=AB4(2,1);
        DHU(3)=(UTG(4)-UTG(1))/h(k);DHV(3)=(VTG(4)-VTG(1))/h(k);
        % Buoc cuoi
        VPC= h(k)*funt2(ttg(4),UTG(4),VTG(4))/(6*UTG(4))+ UTG(1)+ ttg(4)*VTG(1) + (h(k)/6)*(DHU(1)+ 2*DHU(2) + 2*DHU(3) + ttg(4)*(DHV(1)+ 2*DHV(2)+ 2*DHV(3))); 
        MTC(1,1)=1; MTC(2,2)=-1;MTC(1,2)=t(i+1);MTC(2,1)= exp(-ttg(4));
        ABC=[UTG(4);VTG(4)]- inv(MTC)*(MTC*[UTG(4);VTG(4)]-[VPC;VP4]);
        U(i+1)=ABC(1,1);V(i+1)=ABC(2,1);
        UCX(i+1)=exp(t(i+1));    VCX(i+1)=sin(t(i+1));
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
%ketthuc = clock;
%ketthuc - batdau
function D=funt2(t,U,V)
D=U*V*exp(t) + exp(2*t)+ t*cos(t)*exp(t) - sin(t)*exp(2*t);
