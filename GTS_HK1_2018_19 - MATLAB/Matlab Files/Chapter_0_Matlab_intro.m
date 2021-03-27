% Lecture 1. Simple Matlab command 
% Copyright: Phi Ha (18/03/2020)


clear all; close all; clc
% quit; exit; help; ans; edit; 

%% Log base
a = 3
b = 2
a - b
a + b;
a * b
a.^b

e = exp(1)
log2(2)
log(e)
log10(10000)

%% Loops and condition 
 for i=1:5:60
     i
     if (i>30) && (i<50)
       disp('Thay Phi dep trai vai chuong')
     elseif i<= 30
       disp('Lop minh that nhieu giai dep')
     else
       disp('Lop minh that lam gai xinh')
     end
 end    
     
%%
i = 100;
while i>50
    if mod(i,3) == 1
        i
    elseif mod(i,3) == 2
        disp('I love you')
    else 
        disp('Wanna kiss me?')
    end
    i = i - 1;
end   

%% Define function - inline function - anonymous function and function handle 
y = @(x) x .^3;
z = @(x) x * x + 1 ;
y = inline('x.^3')

y(3)
z(4)
feval(y,3)

x = 1:1:4;
[y,z] = luy_thua(x)

%% Simple matrix operation
A = [1 2 3; 3 4 5; 6 7 8]
A(1,1)

disp('Co cua ma tran A la ')
size(A)

disp('Dinh thuc cua A la ')
det(A)

disp('hang cua A la ')
rank(A)

B = eye(3)
C = zeros(3)
D = ones(3)

e = ones(1,4)
E = diag(e)
diag(E)
%% Random number and matrix for test 
x = rand
A = rand(3,3)


%% Plot of function, multifunction, multiplots, hold on hold off
x = 1:1:10;
y = x.^2;

figure(1); clf;
plot(x,y,'b-','linewidth',5)
title('Do thi ham y = x^2')
xlabel('x')
ylabel('x^2')
legend('x^2')

figure(2); clf;
fplot(@(x)sin(x),[-10 10],'rs')
legend('sin(x)')
grid on

figure(3); clf;
subplot(2,2,1)
plot(x,y,'b-')

subplot(2,2,4)
fplot(@(x)sin(x),[-10 10],'rs'); hold on
fplot(@(x)cos(x),[-10 10],'b-.*'); hold off
legend('sin(x)','cos(x)')
grid on

 

