# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 10:21:27 2020

Solving Linear algebraic systems

1. Gaussian elimination, LU, PLU, Cholesky decompositions.
2. Iterative methods: Jacobi, Gauss-Seidel
 
@author: phiha
"""

# condtion number of the matrix A
def kA(A):   
    ivA = la.inv(A)  # Find the inverse of a matrix A
    value = la.norm(A,2) * la.norm(ivA,2)    
    return value             

import numpy as np    
import numpy.linalg as la  # import the module linear algebra in numpy
from scipy.linalg import lu, cholesky, svd,qr

#A = np.array([[1, 0],[0, 1]])
#A = np.eye(10)

A = np.random.rand(3,3)
kA(A)

from scipy.linalg import tril
print(tril(A))

#%%

import scipy.linalg as sla
A = sla.hilbert(10)
P,L,U = lu(A)   
print('P is: ',P)
print('L is: ',L)
print('U is: ',U)    

#import time
#print("Pause for some sec")
#time.sleep(5.5)    # pause 5.5 seconds

A = sla.hilbert(10)
print('Phân tích Cholesky trong scipy.linalg cho 1 ma trận tam giác trên U1 thỏa mãn U1.T * U1 = A')
U1 = cholesky(A)
#print('U1 is: ',U1)
#print('Chuyển vị của U1 là: U1.T = ',U1.T)


print('Phan tich SVD cua A')
U,S,V = svd(A)
#print('S is: ',S)
#print('U is: ',U)


#%%   QR  - infact scipy and numpy give different Q and R, why?

print('===================================================================')
A = np.array([[1,1,1,1,1],[-2,-1,0,1,3],[4,1,0,1,9]]).T
print(A)
b = np.array([-1,0,2,3,4])
print(b)


from numpy.linalg import qr

print('Phan tich QR cua A')
Q,R = qr(A)
print('Q is: ',Q)
print('R is: ',R)
R1 = R[0:3,0:3]

from numpy.linalg import solve
from numpy import dot

a = solve(R1,dot(Q.T,b))
print(a)

Q1 = np.zeros((3,6));
Q2 = np.random.rand(3,6) ;
A = np.concatenate((Q1,Q2))
U,S,V = svd(A)
#print('S is: ',S)

#%%

print('===================================================================')

A = np.array([[2.34, -4.10, 1.78],[1.98, 3.47, -2.22],[2.36, -15.17, 6.81]])
print(A)
b = np.array([0.02, -0.73, -6.63])
print(b)

P,L,U = lu(A)   
print('P is: ',P)
print('L is: ',L)
print('U is: ',U)    





