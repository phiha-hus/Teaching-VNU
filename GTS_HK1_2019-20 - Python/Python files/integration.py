# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 13:24:10 2020

@author: phiha
"""

#%% 
"""
Phuong phap hinh thang composite
de xap xi tich phan ham f tren doan [a,b]
"""
import numpy as np

def trap_int(f,a,b,n):        
    h = (b-a)/n;
    x = np.linspace(a,b,n);    
    Sn = h/2 * ( 2 * sum(f(x)) - f(a) - f(b) );    
    return Sn

def trap(f,a,b,tol):
    n = 2;
    Sn = trap_int(f,a,b,n)
    S2n = trap_int(f,a,b,2*n)
    
    while abs(S2n-Sn)>tol:
        n *= 2;
        Sn = S2n;
        S2n = trap_int(f,a,b,2*n)
    return S2n, n

# Test composite trapezoidal rule
def f(x): return 10*x**9
(a,b,tol) = (0,1,1e-9);

err_0 = abs(trap_int(f,a,b,2)-1)

for i in range(2,13):
    n = 2**i
    Sn = trap_int(f,a,b,n)
    err_n = abs(Sn-1)
    ratio = err_0/err_n
    err_0 = err_n
    print(ratio)


#%% 
import numpy as np
def trap2(f,a,b,tol):
    n = 2;
    Sn = trap_int(f,a,b,n)
    h = (b-a)/n
    x = np.linspace(a,b,2*n)
    x2 = x[1:2*n:2]
    S2n = Sn/2 + h/2 * sum(f(x2)) ; 
    
    while abs(S2n-Sn)>tol:
        n *= 2;
        Sn = S2n;        
        h = (b-a)/n
        x = np.linspace(a,b,2*n)
        x2 = x[1:2*n:2]
        S2n = Sn/2 + h/2 * sum(f(x2)) ; 
    return S2n, n
    
# Test composite trapezoidal rule
def f(x): return 100*x**99    
(a,b,tol) = (0,1,1e-6);
I, n = trap2(f,a,b,tol)  
print(I)
print(n)  

#%% 
"""
Phuong phap Simpson composite 
de xap xi tich phan ham f tren doan [a,b]
"""
import numpy as np
#import numpy.linalg as la

def simpson_int(f,a,b,n):    
    if n%2 == 1:
        n += 1
        
    h = (b-a)/n;    
    x = np.linspace(a,b,n+1);    
    
    # Notice start:stop:step
    x1 = x[1:n:2]; 
    #print(x1)
    
    x2 = x[0:n+1:2];    
    #print(x2)
    
    I = - f(a) + 4 * sum(f(x1)) + 2 * sum(f(x2)) - f(b) ;
    I *= h/3;
    
    return I

def tichphan_simpson(f,a,b,tol):
    n = 2;
    I_old = simpson_int(f,a,b,n);
    #print(I_old)
    
    n *= 2;
    I_new = simpson_int(f,a,b,n);
    #print(I_new)
    
    while abs(I_new-I_old)>tol:        
        I_old = I_new;        
        
        n *= 2;        
        I_new = simpson_int(f,a,b,n);        
    return I_new, n
 
# Test Simpson
def f(x): return 100*x**99
I = simpson_int(f,0,1,2)  
print(I)

(a,b,tol) = (0,1,1e-6);
I, n = tichphan_simpson(f,a,b,tol)  
print(I)
print(n)

