% Test Example em Linh nho
% Symbolic matrix - Matlab R2012b
% Page 20
% Notation: lb: lambda, si: sigma
% The matrix will be divided into 4 blocks for easier implementation
% S = [A B;C D]

clear all; close all; clc
syms lb lb1 lb2 lb3 si1 si2 si3 x1 x2 y1 y2 z1 z2

% Block A
A = sym(zeros(6,7));
for i=1:4
    A(i,i) = 1;
    A(i+2,i) = -1;
end
A(1,5) = lb1 * x1;
A(2,5) = lb1 * x2;
A(3,6) = lb2 * y1;
A(4,6) = lb2 * y2;
A(5,7) = lb3 * z1;
A(6,7) = lb3 * z2;
A

% Block B
B = sym(zeros(6,7));
for i=1:4
    B(i,i) = 1;
    B(i+2,i) = -1;
end
B(1,5) = j * lb1 * x1;
B(2,5) = j * lb1 * x2;
B(3,6) = j * lb2 * y1;
B(4,6) = j * lb2 * y2;
B(5,7) = j * lb3 * z1;
B(6,7) = j * lb3 * z2;
B

% Block C
C = sym(zeros(6,5));
C(1,1) = si1 * lb1; 
C(2,2) = si1 * lb1;
C(3,1) = si2 * lb2; 
C(4,2) = si2 * lb2;
C(5,1) = si3 * lb3; 
C(6,2) = si3 * lb3;

C(1,3) = si1 * lb1 * lb1 * x2; C(1,4)=x2/lb1;
C(2,3) = -si1 * lb1 * lb1 * x1; C(2,4)=-x1/lb1;
C(3,3) = si2 * lb2 * lb2 * y2; C(3,4)= -y2/lb2; C(3,5)=y2/lb2;
C(4,3) = -si2 * lb2 * lb2 * y1; C(4,4)=y1/lb2; C(4,5)=-y1/lb2;
C(5,3) = si3 * lb3 * lb3 * z2; C(5,5)=-z2/lb3;
C(6,3) = -si3 * lb3 * lb3 * z1; C(6,5)=z1/lb3;
C

% BloDk D
D = sym(zeros(6,5));
D(1,1) = j * si1 * lb1; 
D(2,2) = j * si1 * lb1;
D(3,1) = j * si2 * lb2; 
D(4,2) = j * si2 * lb2;
D(5,1) = j * si3 * lb3; 
D(6,2) = j * si3 * lb3;

D(1,3) = -si1 * lb1 * lb1 * x2; D(1,4)=-j * x2/lb1;
D(2,3) = si1 * lb1 * lb1 * x1; D(2,4)=j * x1/lb1;
D(3,3) = -si2 * lb2 * lb2 * y2; D(3,4)= j * y2/lb2; D(3,5)=-j * y2/lb2;
D(4,3) = si2 * lb2 * lb2 * y1; D(4,4)=-j * y1/lb2; D(4,5)=j * y1/lb2;
D(5,3) = -si3 * lb3 * lb3 * z2; D(5,5)= j * z2/lb3;
D(6,3) = si3 * lb3 * lb3 * z1; D(6,5)= -j * z1/lb3;
D

S= [A C;B D]

det(S)

simplify(ans)
