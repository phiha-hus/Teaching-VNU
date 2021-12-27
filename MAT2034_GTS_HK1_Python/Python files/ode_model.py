# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 13:36:32 2020

@author: phiha
"""
import numpy as np
from ivp_ode import trap_imp

(T0,I0,V0) = (1000,10,20) ; 
x0 = np.array([T0, I0, V0]).T
(beta,delta,c,p) = (0.3,0.2,0.3,0.1);
t0 = 0; tf = 10;


def F(t,X):    
    value = np.array([-beta * X[0] * X[2], beta * X[0] * X[2] - delta * X[1] , 
                      p * X[1] - c * X[2] ]).T ;   
    return value

n = 1000;
I, x = trap_imp(F,t0,tf,x0,n) ;
x = np.array([x])

import matplotlib.pyplot as plt
plt.plot(I,x[:,:,0].T,I,x[:,:,1].T,I,x[:,:,2].T) ;
plt.legend(['T','I','V'])
plt.grid() ; 
plt.show()

