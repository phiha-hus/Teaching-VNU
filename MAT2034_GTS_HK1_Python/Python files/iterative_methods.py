# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 14:18:23 2020

@author: phiha
"""


from numpy import array, dot
import scipy.linalg as sla

A = array([[3,2, 1],[1, -5, 3],[1, 2, -8]])

b = array([6,-1,-5])

#P =  array([[3, 0, 0],[0, -5, 0],[0, 0, -8]])  # Jacobi
P =  array([[3, 0, 0],[1, -5, 0],[1, 2, -8]]) # Gauss-Seidel
N = A - P; tol = 1e-9;

#x0 = array([1.1,0.9,0.8]); 
x0 = array([0,0,0]); 
xold = x0;
nmax = 1000;

for i in range(1,nmax):
    xnew = sla.solve(P, -dot(N,xold)+b)
    
    
    if (sla.norm(dot(A,xnew)-b,2)<tol):
        print('No by Jacobi method is ',xnew)
        print('i = ',i+1)
        break
    else:
        xold = xnew
    
    