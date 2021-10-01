# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 15:50:41 2019
Algorithm Gauss 1
@author: Đỗ Văn Chung
"""
def Gauss_1(f,a,b,k):
    from numpy import array, sqrt
    
    if (k <= 3):
        n = 2;
        t = array([1/sqrt(3), -1/sqrt(3)])
        w = array([1, 1])            
    elif (4<=k) and (k<=5):
        n = 3;
        t = array([0, sqrt(0.6), -sqrt(0.6)])       
        w = array([8/9, 5/9, 5/9])
    
    x = (b+a)/2 + (b-a)/2 * t
    temp = w * f(x)
    I = (b-a)/2 * sum(temp)
    
    return I

# Test
def f(x): return x
a = 0; b = 1;
k = 3    
I = Gauss_1(f,a,b,k)
print(I)