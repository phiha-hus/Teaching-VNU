# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 06:52:15 2019
module Solving Algebraic Equation
function included: tim_khoang_nghiem, phan_doi(bisection), 
lap_don (fixed point iteration), newton
@author: Phi Ha
"""

def tim_khoang_nghiem(f,a,b,dx):
    from numpy import sign
    #import math
    #import error
    
    x1 = a; x2 = a+dx; 
    f1 = f(x1); f2 = f(x2);
    
    if sign(f1) == sign(f2) and x2<b :
        #error.err("khong phai khoang nghiem")        
        x1 = x2; x2 = x1+dx;
        f1 = f2; f2 = f(x2);
    elif x2>b or x2==b:
        print("Co the khong co nghiem trong khoang dang xet")
    else:
        print("Da tim thay khoang nghiem")
        exit
    
    return x1, x2
        
def f(x):
    return x**3
    
#a,b,dx = (1,2,1e-6)
#x1, x2 = tim_khoang_nghiem(f,a,b,dx)
#print(x1,x2)

#########################################################################
# Phan 2: Phuong phap phan doi

def bisection(f,a,b,tol):    
    from numpy import sign
    from math import log2
    
    x1 = a; x2 = b;
    f1 = f(a); f2 = f(b);
    
    imax = log2((b-a)/tol);
    print("imax la --- ",imax)
    i = 1 ;
    
    while i < imax :
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

def f(x):
    from math import sin
    return sin(x)-x - 0.01

#a,b,tol = (-0.477777,0.00001,1e-12)
#c, i = bisection(f,a,b,tol)
#print(c)
#print(i)

#########################################################################
# Phan 3: ve hinh 
def f2(x):
    return x**3 - 6 * x**2 + 11*x - 6

def f3(x):    
    return x**2-2*x+1

import matplotlib.pyplot as plt
import numpy as np
#from math import sin

I = np.linspace(0,5,100);

#plt.subplot(1,2,1)
plt.plot(I,f2(I),'-b',I,f3(I),'-.r',linewidth=3)
plt.grid()
plt.title("Do thi ham so f(x)=x**3 - 6 * x**2 + 11*x - 6")
plt.xlabel("x")
plt.ylabel("y")
plt.legend(("f2(x)","f3(x)"))
plt.axis((0,10,-0.5,0.5))
plt.show(1)

def f4(x):
    from scipy import tan
    return tan(x)-x

#I = np.linspace(0,100,1000);
#plt.plot(I,f4(I),linewidth=3)
#plt.grid()
#plt.axis((0,10,-0.5,0.5))


#for i_x, i_y in zip(I,f4(I)):
#    plt.text(i_x, i_y, '({}, {})'.format(i_x, i_y))

#plt.show()

            
#########################################################################
# Phan 4: Phuong phap lap don
# Giai phuong trinh f(x)=0 nhung duoc viet duoi dang moi  x= g(x)
"""
g : ham so, x0: gia tri ban dau, bat ky
nmax : so luong toi da cac phep lap, tol: sai so
x0 : nghiem, i: so luong phep lap thuc te de dat duoc sai so tol
"""
def lapdon(g,x0,nmax,tol):    
    i = 1;
    
    for i in range(nmax):
        x1 = g(x0); 
        g1 = g(x1);
        if abs(x1-g1)>tol :
            x0 = x1;
            x1 = g(x1);
            g1 = g(x1);
            i += 1;
        else:
            print("No la ---",x1);
            break
        if i>nmax:
            print("Co the phep lap don nay 0 cho ket qua")
            x1 = None;
    return x1, i

# Test code
from numpy import exp, arctan

#def g(x): return exp(-x)  # Bai tap 6a Exercise sheet No 34
def g(x): return 1+arctan(x)  # Bai tap 6b Exercise sheet No 34

#x0,nmax,tol = (0,100,1e-6)
#c, i = lapdon(g,x0,nmax,tol)
#print("No la ---",c)
#print("So lan lap la ---",i)


#########################################################################
"""
Phan 5: Phuong phap Newton giai phuong trinh vo huong
Giai phuong trinh vo huong f(x)=0 su dung phuong phap Newton
f : ham so, df: dao ham, x0: gia tri ban dau, bat ky (co the tim duoc bang 
phuong phap phan doi)
nmax : so luong toi da cac phep lap, tol: sai so
x0 : nghiem, n: so luong phep lap thuc te de dat duoc sai so tol
Chu y: x0 la 1 so chu 0 phai vector 
"""

def newton(f,df,x0,nmax,tol):
    x1 = x0
    for n in range(0,nmax):
        fx1 = f(x1)
        if abs(fx1) < tol:
            print('Tim duoc nghiem sau',n,'lan lap.')
            return x1
        Df = df(x1)
        if Df == 0:
            print('Khong tim thay nghiem vi dao ham = 0')
            return None
        x1 = x1 - fx1/Df                
    print('Vuot qua so lan lap, khong tim thay nghiem.')
    return None

# Test code
#x0,nmax,tol = (1,10,1e-10)
#f = lambda x: x**3 - x**2 - 1
#df = lambda x: 3*x**2 - 2*x
#approx = newton(f,df,x0,nmax,tol)
#print(approx)
#
#
#f = lambda x: x**(1/3)
#df = lambda x: (1/3)*x**(-2/3)
#approx = newton(f,df,0.1,100,1e-2)
#print(approx)

#########################################################################
"""
Phan 6: Phuong phap Newton giai he phuong trinh
Giai phuong trinh vo huong f(x)=0 su dung phuong phap Newton
f : ham so, df: dao ham, x0: gia tri ban dau, bat ky (co the tim duoc bang 
phuong phap phan doi)
nmax : so luong toi da cac phep lap, tol: sai so
x0 : nghiem, n: so luong phep lap thuc te de dat duoc sai so tol
Source: HPL
"""
import numpy as np

def Newton_system(F, J, x, nmax, tol):
    """
    Solve nonlinear system F=0 by Newton's method.
    J is the Jacobian of F. Both F and J must be functions of x.
    At input, x holds the start value. The iteration continues
    until ||F|| < tol.
    """
    F_value = F(x)
    F_norm = np.linalg.norm(F_value, ord=2)  # l2 norm of vector
    iteration_counter = 0
    while abs(F_norm) > tol and iteration_counter < nmax:
        delta = np.linalg.solve(J(x), -F_value)
        x = x + delta
        F_value = F(x)
        F_norm = np.linalg.norm(F_value, ord=2)
        iteration_counter += 1

    # Here, either a solution is found, or too many iterations
    if abs(F_norm) > tol:
        iteration_counter = -1
    return x, iteration_counter



#from numpy import cos, sin, pi
#
#def F(x):
#    return np.array([x[0]**2 - x[1] + x[0]*cos(pi*x[0]),
#                     x[0]*x[1] + exp(-x[1]) - x[0]**(-1)])
#
#def J(x):
#    return np.array([[2*x[0] + cos(pi*x[0]) - pi*x[0]*sin(pi*x[0]), -1],
#                      [x[1] + x[0]**(-2), x[0] - exp(-x[1])]])
#
#expected = np.array([1, 0])
#tol = 1e-4
#x, n = Newton_system(F, J, x=np.array([2.0, -1.0]), nmax = 100, tol=0.0001)
#print(n)
#print(x)
#error_norm = np.linalg.norm(expected - x, ord=2)
#assert error_norm < tol, 'norm of error =%g' % error_norm
#print('norm of error =%g' % error_norm)

#############################################################################
def main():
    """
    Module of chapter 3 - solving nonlinear equations and systems
    """
    if __name__ == "__main__":
        main()