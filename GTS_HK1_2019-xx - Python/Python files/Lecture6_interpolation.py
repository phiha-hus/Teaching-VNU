# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 13:58:33 2020
Interpolation in Python

@author: phiha
"""

import scipy as sp
import numpy as np
from scipy.linalg import solve
from numpy import power

def lagrange_interp(xData,yData):
    
    m = xData.shape[0];
    
    A = np.zeros((m,m));
    
    A[:,0] = np.ones(m)
    A[:,1] = xData ; 
    
    for i in range(2,m):
        A[:,i] = power(xData,i)
    
    p = solve(A,yData)
    
    return p

# test
x = np.array([0, 1, 2, 3])
print(x.T)
y = np.array([1, 5, 10, 17])
p = lagrange_interp(x,y)
print(p)
z = sp.polyval(p,0.5)
print(z)


#%%
"""
Day la noi suy tuyen tinh tung khuc - 0 phai noi suy Lagrange hay Newton
"""
x = np.array([0, 1, 2, 3])
y = np.array([1, 5, 10, 17])
xeval = np.array([0.5])
yeval = sp.interp(xeval,x,y)
print(yeval)



