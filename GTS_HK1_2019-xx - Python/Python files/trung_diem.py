# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 15:39:09 2019

@author: Đỗ Văn Chung
"""

def trung_diem(f,a,b,tol):
    #from numpy import linspace
    n = 10;
    from tinh_tich_phan import trung_diem_1
    In = trung_diem_1(f,a,b,n);
    
    n_new = n * 2;
    In_new = trung_diem_1(f,a,b,n_new);
    
    while (abs(In_new-In)>tol):
        n = n_new;
        In = In_new;
        
        n_new = n * 2;
        In_new = trung_diem_1(f,a,b,n_new);
    
    print(n)
    print(In)
 
def f(x): return x**3
a = 0; b = 1; tol = 1e-6
trung_diem(f,a,b,tol)
       
    