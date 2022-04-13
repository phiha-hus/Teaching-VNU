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
%% Matrix commands
A = [1 2;3 4]
A*A
A.^2
power(A,2)
A**2 % Error

%% Logic
2 == 3
2 ~= 3
2 >= 3
2 <= 4


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

%% Plot in 2D of function, multifunction, multiplots, hold on hold off
clear all; close all; clc

%x = -10 : 0.1: 10;
%y = x.^2 - 4*x + 4;

X = linspace(-10,10,100);
y = @(x) x.^2 - 4*x + 4;
Y = y(X);
Z = 1 - X.^2;

% figure(10); clf;
% xlabel('Nhan truc hoanh')
% ylabel('Nhan truc tung')
% title('Ve do thi cua 2 ham so ..... ')
% plot(X,Y,'r -',X,Z,'b s')
% legend('y','z')
% grid on
% 
% %
% figure(1); clf;
% plot(X,Y,'b-','linewidth',5)
% title('Do thi ham y = x^2')
% xlabel('x')
% ylabel('x^2')
% legend('x^2')

%%
% figure(2); clf;
% fplot(@(x)sin(x),[-10 10],'rs')
% legend('sin(x)')
% grid on
% 
% figure(3); clf;
% subplot(2,2,1)
% plot(X,Y,'b-')
% 
% subplot(2,2,4)
% fplot(@(x)sin(x),[-10 10],'rs'); hold on
% fplot(@(x)cos(x),[-10 10],'b-.*'); hold off
% legend('sin(x)','cos(x)')
% grid on

figure(4); clf;
subplot(1,2,1)
xlim([-1,1])
ylim([-1,1])
ezplot('x^2 + y^2 - 4 = 0')
grid on

subplot(122)
% x = @(t) sin(3*t) * cos(t)/(t+pi)
% y = @(t) sin(3*t) * sin(t)/(t+pi)
% t = [0, 4 *pi]
ezplot('sin(3*t)* cos(t)/(t+pi)','sin(3*t)* sin(t)/(t+pi)',[0, 4 *pi])
grid on

%% Plot in 3D
clear all; close all; clc

figure(5); clf;
subplot(121)
T = linspace(0,100,1000);
plot3(sin(T),cos(T),T,'b -')

subplot(122)
x = 0:0.05:2;
y = -1:0.2:1;
[X,Y] = meshgrid(x,y);
Z =  (Y+1).*cos(2*pi*X.^2) + (X+1).*cos(2*pi*Y.^2);
%surf(X,Y,Z)
mesh(X,Y,Z)
xlabel('x'); ylabel('y'); zlabel('z')

title('Surface view')

print -f -dpng 'Figure 5'








