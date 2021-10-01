# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 16:02:37 2019

@author: Đỗ Văn Chung
"""

from numpy.linalg import qr, solve
from numpy import array, dot, shape,ones,power,append

T = array([1,2,3,4,5,6,7])
A1 = append(ones((1,7)), [T, power(T,2)],axis=0)
A = A1.T

V1 = array([2.31, 2.01, 1.8, 1.66, 1.55, 1.47, 1.41])
V = V1.T

print("Cach 1: phuong trinh chinh tac")
x = solve(dot(A.T,A),dot(A.T,V))
print(x)

print("Cach 2: QR")
Q,R = qr(A)
print(Q)
print(R)
x = solve(R,dot(Q.T,V))
print(x)
