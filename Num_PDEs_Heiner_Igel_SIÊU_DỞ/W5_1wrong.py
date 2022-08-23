# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 14:44:48 2021

@author: phiha
"""

from numpy import linspace
from numpy import exp, pi, arange

# [...]
# Basic parameters
nx = 128
x0 = pi
def fourier_derivative(f, dx):
    # Length of vector f
    nx = f.size
    # Initialize k vector up to Nyquist wavenumber
    kmax = pi/dx
    dk = kmax/(nx/2)
    k = arange(float(nx))
    k[:nx/2] = k[:nx/2] * dk
    k[nx/2:] = k[:nx/2] - kmax
    # Fourier derivative
    ff = 1j*k*fft(f)
    df = ifft(ff).real
    return df
# [...]
# Main program
# Initialize space and Gauss function (also return dx)
x = linspace(2*pi/nx, 2*pi, nx)
dx = x[1]-x[0]
sigma = 0.5
f = exp(-1/sigma**2 * (x - x0)**2)
# Calculate derivative of vector f
df = fourier_derivative(f, dx)
# [...]