
# -*- coding: utf-8 -*-
"""
Code for polynomial division
"""
def tim_thuong(a,r):
# ham thuc hien phep chia da thuc, trong do
    # a la mang luu cac he so cua p(x)
    # r la so ma (x+r) la so chia
    n = len(a);
    b = np.zeros(n-1);
    b[0] = a[0];
    
    for i in range(1,n-1):
        b[i] = a[i] - r * b[i-1]
    
    s = a[n-1] - r * b[n-2];
    
    print('He so cua thuong la ---',b)    
    print('So du la --- ',s)


    return b, s
       
# test ham tim_thuong
import numpy as np
A = np.array([1, 3, 3]);
r = 1;
b, s = tim_thuong(A,r)
print(b)

c = np.ones(5)
print(c)


def giai_tamthuc(A):
    n = len(A)
    if n!=3:
        print("Loi 0 phai tam thuc")
    elif A[0]==0:
        print("Khong phai tam thuc")
    else:
        delta = A[1]**2 - 4 * A[0]*A[2];
        if (delta<0):
            print("Vo nghiem")
        elif (delta==0):
            x = A[1]/2;
            print(x)
        else:
            from math import sqrt
            x1 = (-A[1] - sqrt(delta))/(2*A[0]);
            x2 = (-A[1] + sqrt(delta))/(2*A[0]);
            print(x1)
            print(x2)

# test
import numpy as np
A = np.array([1,3,2])
giai_tamthuc(A)
