# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 16:51:05 2019
@author: Phi Ha
"""

def bt8(x0,nmax,tol):
    from scipy import exp, log
    def f(x): return exp(x)-log(x+2)-2;
    def df(x): return exp(x) - 1/(x+2)
    
    for n in range(nmax):
        x1 = x0 - f(x0)/df(x0)
        if abs(x1-x0)<=tol:
            break
        else:
            x0 = x1;
    print(f(x1))    
    return x1, n

# Test
(x0,nmax,tol) = (0,100,1e-5)
x0, n = bt8(x0,nmax,tol)

print(x0)
print(n)
print("Ket thuc bai tap 8")

"""
Bai tap 10 sach cua thay Pham Ky Anh, tim 2 nghiem cua phuong trinh f(x)=0 nam trong lan can cua 2 
Y tuong: 1. Su dung thuat toan tim khoang nghiem de dinh vi 2 khoang nho chua 2 nghiem can tim.
		 2. Su dung Newton voi moi khoang nho do de tim nghiem can tim.
"""
N = 1000;
from numpy import linspace
from scipy import sign

I = linspace(1,3,N);
S = [];

def f(x): return x**4-5*x**3-12*x**2+76*x-79
def df(x): return 4*x**3-15*x**2-24*x+76

for i in range(len(I)-1):
    if sign(f(I[i])) != sign(f(I[i+1])):
        S.append(I[i])
print(S)

(nmax,tol) = (100,1e-6)

for i in range(2):
    x0 = S[i];
    for n in range(nmax):
        x1 = x0 - f(x0)/df(x0)
        if abs(x1-x0)<=tol:
            break
        else:
            x0 = x1;
    print("Nghiem thu ",i+1,"cua chung ta la ",x1)    

