% TITLE: Plot Taylor polynomials for exp(x) about x = 0
%
% This plots several Taylor polynomials and their errors
% for increasing degrees.  The particular function being
% approximated is exp(x) on [-b,b].  First we plot the 
% Taylor polynomials and then we plot the errors.
%
% Initialize
b = input('Give the number b defining the interval [-b,b] ');
h = b/100;
x = -b:h:b;
max_deg = 20;

% Produce the Taylor coefficients for the function exp(x) when
% expanded about the point a = 0.  The coefficients are stored 
% in the array c, which will have length max_deg+1.
c = ones(max_deg+1,1);
fact = 1;
for i = 1:max_deg
  fact = i*fact;
  c(i+1) = 1/fact;
end

% Calculate the Taylor polynomials
pa = polyeval(x,0,c,1);
p2 = polyeval(x,0,c,3);
pb = polyeval(x,0,c,5);
p4 = polyeval(x,0,c,8);

% Calculate the errors in the Taylor polynomials
true = exp(x);
err1 = true-pa;
err2 = true-p2;
err3 = true-pb;
err4 = true-p4;

% Initialize for plotting the Taylor polynomials
hold off
clf
m = min([min(pa),min(p2),min(pb),min(p4),min(true)]);
if m > 0
    m = .9*m;
else
    m = 1.1*m;
end
M = max([max(pa),max(p2),max(pb),max(p4),max(true)]);
if M < 0
    M = .9*M;
else
    M = 1.1*M;
end
axis([-b,b,m,M])
hold on

% Plot the Taylor polynomials
plot(x,true,x,pa,':',x,p2,'--',x,pb,'-.');
title('Taylor approximations of e^x')
plot([-b,b],[0,0])
plot([0 0 ],[m M/1.2])
axis('equal')
legend('exp(x)','degree = 1','degree = 3','degree = 5',0)
pause

% print -deps exp2d.eps
% print

% Initialize for plotting the errors
hold off
clf
m = min([min(err1),min(err2),min(err3),min(err4)]);
if m > 0
    m = .9*m;
else
    m = 1.1*m;
end
M = max([max(err1),max(err2),max(err3),max(err4)]);
if M < 0
    M = .9*M;
else
    M = 1.1*M;
end
axis([-b,b,m,M])
hold on

% Plot the errors
plot(x,err1,x,err2,':',x,err3,'--',x,err4,'-.');
title('Errors in Taylor approximations of e^x')
legend('degree = 1','degree = 3','degree = 5','degree = 8',0)
hold off

% print -deps experr2d.eps
% pause
% print
