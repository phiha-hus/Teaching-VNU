# -*- coding: utf-8 -*-
"""
Created on Sun Dec 26 00:07:16 2021

@author: DELL
"""

# Test module ivp_ode
import numpy as np   
from ivp_ode import euler_imp, trap_imp, Heun


# Test Verhulst equation
# r, K, N0 = 3, 1000, 100 
# def f(t,N):
#     return r * N * (1- N/K)
# (t0,tf,x0,n) = (0,10,N0,1000)


#g, q = 25, 25
#def f(x,y):
#    return g - q * y

def f(x,y):
    return - 4 * y + x**2

from numpy import exp
def y_exact(x):
    return 31/32 * exp(-4*x) + 1/4 * x**2 - 1/8 * x + 1/32
    

# Change n = 10 and n = 20 to see
(t0,tf,x0,n) = (0,10,1,1000)

# Choose a suitable method
I2, N2 = euler_imp(f,t0,tf,x0,n)
I3, N3 = Heun(f,t0,tf,x0,n)
I4, N4 = trap_imp(f,t0,tf,x0,n)

N2 = abs(np.array([N2]) - y_exact(I2))
N3 = abs(np.array([N3]) - y_exact(I3))
N4 = abs(np.array([N4]) - y_exact(I4))


import matplotlib.pyplot as plt
# plt.plot(I2,N2.T,I3,N3.T,I4,N4.T)
# plt.legend(('euler_imp','Heun','Trap'))
# plt.grid()
# plt.show()

plt.subplot(221)
plt.semilogy(I2,N2.T,'r-')
plt.legend(['euler_imp'])
plt.grid()

plt.subplot(222)
plt.semilogy(I3,N3.T,'b-.')
plt.legend(['Heun'])
plt.grid()

plt.subplot(223)
plt.semilogy(I4,N4.T,'g-')
plt.legend(['Trap'])
plt.grid()

#plt.show()

plt.savefig('Sai_so.png', dpi=150)
plt.savefig('Sai_so.pdf', dpi=150)

