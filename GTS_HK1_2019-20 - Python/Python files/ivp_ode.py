# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 10:19:26 2020

@author: phiha

Module solving differential equations - IVP
"""

def euler_exp(f,t0,tf,x0,n):  
    import numpy as np
    I = np.linspace(t0,tf,n) ;
    h = (tf-t0)/n ;
    x = [x0];    
    xi = x0 ;
    
    for i in range(1,n):        
        x_ip1 = xi + h * f(I[i],xi) ;
        xi = x_ip1 ;          
        x.append(x_ip1)
        #x = np.array([x, xi]);                
    return I, x

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def euler_imp(f,t0,tf,x0,n):  
    import numpy as np
    I = np.linspace(t0,tf,n) ;
    h = (tf-t0)/n ;
    x = [x0];    
    x_old = x0 ;
    
    for i in range(1,n):        
        from scipy.optimize import fsolve
        def f2(xi): return xi - h * f(I[i],xi) - x_old        
        xi = fsolve(f2,x_old)
        #xi = x0
        x_old = xi ;          
        x.append(xi)
        #x = np.array([x, xi]);        
        
    return I, x


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def trap_exp(f,t0,tf,x0,n):
    import numpy as np
    I = np.linspace(t0,tf,n+1) ;
    h = (tf-t0)/n ;
    x = [x0];    
    xi = x0 ;
    
    for i in range(0,n):        
        y_ip1 = xi + h * f(I[i],xi) ;
        x_ip1 = xi + h/2 * ( f(I[i],xi) + f(I[i+1],y_ip1) );
        xi = x_ip1 ;          
        x.append(x_ip1)
        #x = np.array([x, xi]);        
        
    return I, x   
    

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def trap_imp(f,t0,tf,x0,n):
    import numpy as np
    from scipy.optimize import fsolve
    I = np.linspace(t0,tf,n+1) ;
    h = (tf-t0)/n ;
    x = [x0];    
    xi = x0 ;
    
    for i in range(0,n):        
        def f2(x_ip1): return x_ip1 - xi - h/2 * ( f(I[i],xi) + f(I[i+1],x_ip1) )
        x_ip1 = fsolve(f2,xi)
        #xi = x0
        xi = x_ip1 ;          
        x.append(x_ip1)
        #x = np.array([x, xi]);        
        
    return I, x    

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
# Test module ivp_ode
import numpy as np   

# Test Malthus equation
def f(t,x): return t/x

(t0,tf,x0,n) = (0,0.5,2,1000)

# Choose a suitable method
# I, x = euler_exp(f,t0,tf,x0,n)
# I, x = euler_imp(f,t0,tf,x0,n)
I, x = trap_exp(f,t0,tf,x0,n)
#I, x = trap_imp(f,t0,tf,x0,n)

x = np.array([x])
print(x)
#print(type(x))
#print(x.size)

import matplotlib.pyplot as plt
plt.plot(I,x.T)
#plt.plot(I[0:n-1],x.T)
plt.grid()
plt.show()


import matplotlib.pyplot as plt
plt.plot(I,x)
plt.show()
"""
