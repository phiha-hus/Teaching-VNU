# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 14:34:01 2021

@author: phiha
"""

def fourier_derivative(f, dx):
    # Length of vector f
    nx = f.size
    # Initialize k vector up to Nyquist wavenumber
    kmax = pi/dx
    dk = kmax/(nx/2)
    k = arange(float(nx))
    ik = int(nx/2)
    k[:ik] *= dk
    k[ik:] = k[:ik] - kmax
    # Fourier derivative
    from numpy.fft import fft, fftshift
    from numpy.fft import ifft
    
    ff = 1j*k*fft(f)    
    df = ifft(ff).real
    return df


from numpy import linspace
from numpy import exp, pi, arange

# Basic parameters
#nx = 128
nx = 2000
x0 = pi

# Initialize space and Gauss function (also return dx)
x = linspace(2*pi/nx, 2*pi, nx)
dx = x[1] - x[0]
print('dx is: ',dx)

sigma = 0.5
f = exp(-1/sigma**2 * (x - x0)**2)
df = - 2 * (x-x0)**2/sigma**2 * f
#f = (x-pi) **2
#df = 2 * (x-pi)

# Calculate derivative of vector f
df = fourier_derivative(f, dx) - df
print(max(df))