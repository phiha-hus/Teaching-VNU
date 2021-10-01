# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 16:19:43 2019

@author: Đỗ Văn Chung
"""

def Gauss_quad(f,a,b,k,tol):
    nmax = 1e+4
    n = 1;
    
    from numpy import linspace 
    from Gauss1 import Gauss_1
    
    while (n<nmax):
        x = linspace(a,b,n); In = 0;         
        for i in range(n-1):
            In += Gauss_1(f,x[i],x[i+1],k)
        
        x = linspace(a,b,2*n); In_new = 0;
        for i in range(2*n-1):
            In_new += Gauss_1(f,x[i],x[i+1],k)
            
        if abs(In-In_new)<tol:
            break
            return In, n
        else:
            n *= 2;
            
    return In, n

# Test
from numpy import exp    
def f(x): return exp(x)
a = 0; b = 1; tol = 1e-6;
k = 3    
I, n = Gauss_quad(f,a,b,k,tol)
print(I)
print(n)