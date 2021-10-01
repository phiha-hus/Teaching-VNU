% Check how fast the usual multiplication in power function
% in comparison with a modified one.

clear all; close all; clc

%%
result = [];

 %nmin = 3e2; nmax = 4e2;
nmin = 1e3; nmax = 2e3;

for n=nmin:nmax
    x = rand(n,1); p = x; 
    
    tic    
    for i=1:5
      p = p .* p;
    end    
    time1 = toc;
    
    tic
    q = x.^32;
    time2 = toc;
    
    r = x;
    tic
    for i=1:31
        r = x .* r;
    end
    time3 = toc;
    
    result = [result [time2/time1; time3/time1; (time2>time1); (time3>time1)]];    
end

result 

I = linspace(nmin,nmax,nmax-nmin+1);

subplot(1,2,1)
plot(I,result(1,:))
grid on

subplot(1,2,2)
plot(I,result(2,:))
grid on

%%
%n = 100;
%A = zeros(n,n); B=A;
%
%for j = 1 : n 
%  for i = 1 : n 
%     A(i,j)= 1/(i + j-1);
%  end 
%end
%
%for j = 1 : n 
%  for i = 1 : n 
%     B(i,j)= 1.0/real(i + j-1);
%  end 
%end
%
%B-A;

return 