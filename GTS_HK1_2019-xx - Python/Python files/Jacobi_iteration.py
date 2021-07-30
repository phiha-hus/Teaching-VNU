# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 15:24:25 2019

@author: Đỗ Văn Chung
"""

from numpy import array, diag, diagflat,dot
from numpy.linalg import inv, norm

def jacobi(A,b,x0,tol):
    d = diag(A)
    D = diagflat(d)
    R = A-D
    
    G = - dot(inv(D),R)
    f = dot(inv(D),b)
    
    x_old = x0; 
    
    for n in range(100):
        x_new = dot(G,x_old) + f;
        err = norm(dot(A,x_new)-b,1)
        if err<tol:
            print("x can tim la --- ",x_new)
            print("Sai so la ---",err)
            return x_new
            break

    return x_new

##### Test        
A = array([[1,0,0],[0,2,0],[0,0,3]])
b = array([1,2,3])
x0 = array([0,0,0])
tol = 1e-6 

x = jacobi(A,b,x0,tol)
print(x)      
        
    
