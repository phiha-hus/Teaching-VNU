clear all; close all; clc

nmin = 5e1; nmax = 2e2;

T = []; S = linspace(nmin,nmax,nmax-nmin+1);

for n = nmin:nmax
  a = rand(1,n);
  r = rand;
  [y,time1,time2]  = Horner(a,r)
  T = [T [time1; time2; time2/time1; (time2/time1>1)]] 
end

close all; clf;

subplot(2,2,1)
plot(S,T(1,:),'r-.',S,T(2,:),'b*')
title('Compare loop vs. matrix equation in Horner algorithm')
legend('time1','time2')

subplot(2,2,2)
semilogy(S,T(3,:),'r-*')
title('Time2/Time1')

subplot(2,2,3)
plot(S,T(4,:),'b*')