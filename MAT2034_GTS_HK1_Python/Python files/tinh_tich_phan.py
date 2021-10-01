# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 13:29:21 2019

@author: Đỗ Văn Chung
"""

def trung_diem_1(f,a,b,n):
    from numpy import linspace 
    S = linspace(a,b,n+1);  
    #print(S)
    l = len(S)
    #print(l)
    h = (b-a)/n;
    T = S[0:l-1] + h/2  # luoi cac trung diem
    
    # print(T)
    #print(f(T))
    
    #print(sum(f(T)))     
    I = h * sum(f(T))    
    #print(I)
    return I

# Test luon     
def f(x): return x**3

a = 0; b = 1; n = 10;
trung_diem_1(f,a,b,n)

