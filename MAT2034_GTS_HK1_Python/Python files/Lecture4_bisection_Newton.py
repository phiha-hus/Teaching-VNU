# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 08:41:27 2019

@author: MyPC
"""
def bisection(f,a,b,tol):    
    from numpy import sign
    from math import log2
    
    x1 = a; x2 = b;
    f1 = f(a); f2 = f(b);
    
    nmax = log2((b-a)/tol);
    print("nmax la --- ",nmax)
    i =1;
    
    while (i < nmax) :
        if f1 == 0:
            print("No la --- ",x1);
            break
        elif f2 == 0 :            
            print("No la --- ",x2);
            break
        else:
            c = (x1+x2)/2.0; f3 = f(c);
            if sign(f1) != sign(f3):
                x2 = c; f2 = f(c)
            else:
                x1 = c; f1 = f(c)
            i += 1 
    return c, i


def newton(f,df,x0,nmax,epsilon):
    x1 = x0
    for n in range(0,nmax):
        fx1 = f(x1)
        if abs(fx1) < epsilon:
            print('Tim duoc nghiem sau',n,'lan lap.')
            return x1
        Df = df(x1)
        if Df == 0:
            print('Khong tim thay nghiem vi dao ham = 0')
            return None
        x1 = x1 - fx1/Df
        print(x1)                
    print('Vuot qua so lan lap, khong tim thay nghiem.')
    return None
# Test code

x0,nmax,eps = (3,100,1e-6);
f = lambda x: 1-1/x
df = lambda x: 1/x**2

x0, n = bisection(f,1e-3,10,1e-2)
print("x0 tu phan doi la ---",x0)

print("Su dung newton voi x0 tim tu phan doi")
approx = newton(f,df,x0,nmax,eps)
print(approx)
