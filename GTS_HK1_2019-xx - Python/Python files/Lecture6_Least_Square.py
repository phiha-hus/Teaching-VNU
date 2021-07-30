# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 23:10:35 2019
# Least square approximation
# 1st METHOD. Using normal equation
# 2nd METHOD. Using QR decomposition
@author: Phi Ha
"""
import numpy as np
from numpy.linalg import qr,solve
from numpy import array, dot


A = array([[1, 2, 3], [4, 5, 6],[7, 8, 0]])
#A = array([[2., 1., 1.], [1., 3., 2.], [1., 0., 0.]]) 

b = array([1,2,3])

print(A)
print(b)

#%%
"""
1st METHOD. Using normal equation
"""
A1 = dot(A.T,A)
b1 = dot(A.T,b)

print(A1)
print(b1)

x = solve(A1,b1)
print("No la ---", x)

#%%
"""
2nd METHOD. Using QR decomposition
"""

Q,R = qr(A)
print("Q is")
print(Q)

print("R is")
print(R)
R1 = R[:,0:3]

b2 = dot(Q.T,b)
x = solve(R1,b2)
print("No la ---", x)


#%% 
"""
Curve fitting using least square method
Exercise 11 page 143 Kiusalass
"""
import numpy as np
from numpy.linalg import qr,solve
from numpy import array, dot

T = array([79, 190, 357, 524, 690])
k = array([1.0, 0.932, 0.839, 0.759, 0.693])

A = np.ones((3,5))
A[1,:] = T;
A[2,:] = T ** 2
print(A)

#[[1.00000e+00 1.00000e+00 1.00000e+00 1.00000e+00 1.00000e+00]
# [7.90000e+01 1.90000e+02 3.57000e+02 5.24000e+02 6.90000e+02]
# [6.24100e+03 3.61000e+04 1.27449e+05 2.74576e+05 4.76100e+05]]
#
#
#A1 = dot(A,A.T)
#print(A1)
#
#[[5.00000000e+00 1.84000000e+03 9.20466000e+05]
# [1.84000000e+03 9.20466000e+05 5.25238156e+08]
# [9.20466000e+05 5.25238156e+08 3.19648597e+11]]
#
#b1 = dot(A,k)
#print(b1)
#[4.22300000e+00 1.43148900e+03 6.85156395e+05]
#
#x = solve(A1,b1)
#print("No la ---", x)
#[ 1.05258944e+00 -6.80740001e-04  2.30985609e-07]

Q,R = qr(A.T)
print("Q is")
print(Q)

print("R is")
print(R)
R1 = R[:,0:3]

b2 = dot(Q.T,k)
x = solve(R1,b2)
print("No la ---", x)


Q is
[[-0.4472136  -0.58584906  0.51245824]
 [-0.4472136  -0.36083437 -0.13584258]
 [-0.4472136  -0.02229875 -0.56610114]
 [-0.4472136   0.31623686 -0.34143983]
 [-0.4472136   0.65274532  0.53092531]]
R is
[[-2.23606798e+00 -8.22873016e+02 -4.11644909e+05]
 [ 0.00000000e+00  4.93301125e+02  3.78078740e+05]
 [ 0.00000000e+00  0.00000000e+00  8.51676705e+04]]
No la --- [ 1.05258944e+00 -6.80740001e-04  2.30985609e-07]


