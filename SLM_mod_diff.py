# 01/11/2021
# Phase and/or amplitude modulation by SLM/DMD using Fraunhofer diffraction approximation
# Solyanik-Gorgone, creadit to Jonathan George and Paul Hu for partial contributions

from numba import jit
import scipy
import numpy as np
import time
import matplotlib
import matplotlib.pyplot as plt

import pandas
from matplotlib.table import Table

import matplotlib.colors as colors
import math
from math import sqrt
import cmath

padding = 0

complex_type = np.complex128
float_type = np.float64

lambda_0 = 0.5e-6
f = 10e-2
(4*lambda_0*f**3/np.pi)**(1/4)

def propFF(u1,L1,lambda_0,z):
    # Fraunhofer Propagation
    # propagation - Fraunhofer pattern assumes unifom sampling
    # u1 - source plane field
    # L1 - source plane side length
    # lambda_0 wavelngth
    # z - propagation distance
    # L2 - oversvation plane side length
    # u2 oversvation plane field
    # u1 = u1.astype(complex_type)
    #Credit to Jonathan George for this routine re-written in Python
    (M,N)=np.shape(u1)
    L1 = float_type(L1)
    lambda_0 = float_type(lambda_0)
    z = float_type(z)
    M = float_type(M)
    N = float_type(N)
    dx1 = float_type(L1/M)
    k = float_type(2*np.pi/lambda_0)
    L2 = float_type(lambda_0*z/dx1)
    dx2 = float_type(lambda_0*z/L1)
    x2 = np.linspace(-L2/float_type(2),L2/float_type(2)-dx2,int(M))
    (X2,Y2) = np.meshgrid(x2,x2)
    X2 = X2.astype(float_type)
    Y2 = Y2.astype(float_type)
    c = float_type(1)/(complex_type(1j)*lambda_0*z)\
        *np.exp(complex_type(1j)*k/(float_type(2)*z)*(X2**2+Y2**2))
    u2=c*np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(u1)))*(dx1**2)
    return u2

def Phase_OAM(T_0, X, Y, alpha, l, Lambda):
    #Credit to Paul Hu for this function
    r = np.sqrt(X**2 + Y**2)
    theta = np.arctan2(Y, X)
    return (2*(1+np.cos((2.0*np.pi/Lambda) * x-l*theta)))

def screen(amplitude, phase):
    R = np.zeros(np.shape(X)).astype(float)
    R = R.astype(complex_type)

    R.real = np.real(amplitude*np.exp(complex(0, 1) * phase))
    R.imag = np.imag(amplitude*np.exp(complex(0, 1) * phase))

    return R

lambda_0=0.3e-6#wavelength of light
k=2*np.pi/lambda_0
z=1000

#create the grid
dx, dy = 0.05, 0.05
x = np.arange(-144.*.3, 144.*.3, dx)
y = np.arange(-144.*.3, 144.*.3, dy)

X, Y = np.meshgrid(x, y)


#Assign amplitude and/or phase modulation
amplitude = np.ones(np.shape(X)).astype(float)
phase = Phase_OAM(1., X, Y, 1., 10., 1.)

R = screen(amplitude, phase)

#Propagate from the screen in free space
u2_ff = propFF(R,len(R),lambda_0,z)

plt.figure(figsize=(8,7))
original=np.real(u2_ff*np.conj(u2_ff))
plt.pcolormesh(X,Y,original)
plt.colorbar()
plt.title('Fraunhofer Propagation from SLM screen')

plt.show()


