# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 14:20:02 2020
Exercise 12 page 196 Kiusalass
@author: phiha
"""

import numpy as np
import numpy.linalg as la

R = 90 # mm

def x(R,theta):
    from numpy import sin, cos, sqrt
    value = R * (cos(theta) + sqrt(6.25 - sin(theta) * sin(theta)))
    return value

theta = np.linspace(0,180,180)

xData = x(R,theta)

tData = theta/3e+4
h = tData[1] - tData[0]

n = len(xData)

aData = xData[0:n-2] - 2 * xData[1:n-1]+ xData[2:n]
print(aData/h*h)

import matplotlib.pyplot as plt
plt.semilogy(3600 * tData[0:n-2],aData/1000)
plt.show()