# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 10:14:21 2020

@author: phiha
"""

#%%
print("======================================================================")
"""
Example No. 1
"""
import scipy.sparse as sp
from numpy import array
import numpy as np
from numpy.random import rand
from numpy.linalg import solve, norm

m = 2; n = 3;
k = 3
data = np.round(rand(m,n))  # de nhin cho no de thoi ma
print(data)

mtx = sp.lil_matrix((k*n, k*n))
"""
assign the data using fancy indexing: using row based
row index 0:m
column index 1:(n+1)    
"""
mtx[:m, 1:(n+1)] = data   


A1 = mtx.toarray()
A2 = mtx.todense()
print(A1)
print(A2)

#np.reshape(A1,(k*n,k*n))
np.reshape(A1,(9,9))
print('size of A1 is ',np.size(A1))

print(np.size(A2))

b = rand(1,k*n)
print('size of b is ',np.size(b))
"""
error while using np.linalg.solve, I don't know why
"""
# x = solve(A1,b) 

#%%
print("======================================================================")

from numpy.random import rand
from time import time
from numpy.linalg import solve

N = 10000
A = rand(N,N)
b = rand(N)

tic = time()
x = solve(A,b)
toc = time() - tic
print("Processing time for solving system of size %f is %f second"%(N,toc))


#Processing time for solving system of size 10000.000000 is 17.253743 second

#%%
print("======================================================================")
"""
Example No. 1
"""
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from numpy.random import rand
from time import time
import numpy as np

N = 1000000;
A = lil_matrix((N, N),dtype=np.float64)
A[0, :100] = rand(100)
A[1, 100:200] = A[0, :100] * 2

A.setdiag(rand(N))

A = A.tocsr()
b = rand(N)

tic = time()

x = spsolve(A, b)

toc = time() - tic
print("Processing time for solving system of size %f is %f second"%(N,toc))

#Processing time for solving system of size 1000000.000000 is 0.692606 second