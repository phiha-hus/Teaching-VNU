clear all; close all; clc
syms xi1 a1 b1

A = [1 xi1 a1 b1
     1 xi1^2 a1^2 b1^2
     1 xi1^3 a1^3 b1^3
     1 xi1^4 a1^4 b1^4]
 
det(A)
 
simplify(det(A))
 
factor(det(A))
