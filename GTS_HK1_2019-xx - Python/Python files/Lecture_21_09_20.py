# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 10:34:48 2020

@author: phiha
"""

# Test code
from math import sqrt 
from numpy import power, linspace
import matplotlib.pyplot as plt
from time import time
from Lecture3 import newton, bisection

global L, h

L = 16.0; h = 1.5;
x0,nmax,tol = (0.05,100,1e-9)


def f(x):
    global L
    from scipy import cosh
    value = 1/x * (cosh(x*L/2)-1) - h;
    return value

def df(x):
    global L
    from scipy import cosh, sinh
    value = -1/(x**2) * (cosh(x*L/2)-1) + 1/x * sinh(x*L/2) * L/2 ;
    return value

T = linspace(-0.5,0.5,1000)
plt.plot(T,f(T),'-.')
plt.grid()

tic = time()
alpha = newton(f,df,x0,nmax,tol)
#a = -0.2; b = 0.1; tol = 2e-15;
#alpha = bisection(f,a,b,tol)
toc = time() - tic

print('alpha la: ',alpha)
print('sai so cua ham f la: ',abs(f(alpha)))
print('thoi gian tinh toan: ',toc)






##########################################################################################################
def g(x): return 3 * x/4 + 1/power(x,3)

x0,imax,tol = (1.0,100,1e-9)

imax = int(imax)

#from Lecture3 import lapdon
#c, i = lapdon(g,x0,nmax,tol)

# Nghiem chinh xac la 
c = sqrt(2)
i = 0; xi = x0;

while i < imax:
    x_new = g(xi);
    if abs(x_new-c) <= tol:
        imax = i
        break
    else:
        i += 1;
        xi = x_new
        
print(i)

##################################################



