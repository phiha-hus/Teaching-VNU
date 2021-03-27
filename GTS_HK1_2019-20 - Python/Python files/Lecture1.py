# -*- coding: utf-8 -*-
"""
Created on Sun Sep  1 09:41:37 2019
Bai giang so 1
@author: Phi Ha
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from math import sqrt, sin

########################################################################
# Simple way
a = input('Input a: ')
print(a,type(a))
a = eval(a)
print(a,type(a))

# Better way
b = eval(input('Input b: '))



########################################################################
import time
tic = time.time()

for i in range(1,100,1):
    print(sqrt(3*i*i))
    
toc = time.time()
print("Processing time ---- ",toc-tic)


########################################################################
from math import e
print(format(e,'.5f'))
print(format(e,'.5g'))


x = 1./15
sin(x)
print(format(x,'.10e'))
print(format(sin(x),'.10e'))
print(format(x-sin(x),'.10e'))


########################################################################
A = np.array((3,3))
A = np.array([[1,2,3],[4,5,6],[7,8,9]])
print(A)

########################################################################

from math import factorial, exp
n = 1; a = 1;
e = exp(1)


while (abs(e-a)>1e-12):
    a += 1.0/factorial(n)
    print(a)
    
    n += 1
    print(n)


def sq(t):
    return t**2 

import numpy as np
import matplotlib.pyplot as plt

I = np.linspace(-10,10,101);

J = sq(I)

plt.subplot(2,2,1)
plt.plot(I,J,'ro',I,I**2-50,'b-.');
plt.title("Do thi cua ham bac 2")
plt.legend(['x^2','x^3'])

plt.grid()    # khong co cung 0 sao
#plt.show(0)


plt.subplot(2,2,4)
plt.plot(I,J,'ro',I,I**2-50,'b-.');
plt.title("Do thi cua ham bac 2")
plt.legend(['x^2','x^3'])
