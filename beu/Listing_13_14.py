# Reflexion/transmission of quantum wave packet using tridiagonal solver
from math import pi, sqrt, fabs
import numpy as np
import cmath

def Pot(x, a, V0):
    return V0 if fabs(x) <= 0.5e0*a else 0e0

# Initial Gaussian wave packet
def Init(Psi, x, nx, x0, sig, k0):
    # Psi(x,0) = 1/sqrt(sqrt(2*pi)*sig) * exp[-(x-x0)^2/(4*sig^2)] * exp(ikx)
    # x0 - position of center, sig - half-width, k0 - average wave number
    a = 1e0/sqrt(sqrt(2*pi)*sig)
    b =-1e0/(4*sig*sig)
    for i in range(1,nx+1):
        dx = x[i] - x0
        f = a * np.exp(b*dx*dx)
        if (f < 1e-10): f = 0e0
        Psi[i] = f*cmath.exp(-1j*k0*x[i]) # probably np.math is also can be used


# Calculates the probability density Psi2[] of the wave function Psi[]
def ProbDens(Psi, Psi2, nx, hx):
    for i in range(1,nx+1):
        Psi2[i] = abs(Psi[i])*abs(Psi[i])
        # unnormalized probability density
        if (Psi2[i] <= 1e-10): Psi2[i] = 0.0

    # integral by trapezoidal rule    
    PsiNorm = 0.5e0*(Psi2[1] + Psi2[nx])
    for i in range(2,nx):
        PsiNorm += Psi2[i]
    PsiNorm *= hx

    for i in range(1,nx+1):
        Psi2[i] /= PsiNorm

    return PsiNorm # normalized prob. density

# main
a = 5.0     # width of potential barrier
V0 = 50.0   # height of potential barrier
x0 = -20.0  # initial position of wave packet
sig = 1.0    # half-width of packet
k0 =  10.0   # average wave number of packet
xmax = 100.0
hx = 5e-2 # spatial step size
tmax = 5.0
ht = 5e-3 # time step
nout = 40 # output every nout steps

Nx = 2*int((xmax/hx + 0.5) + 1) # odd number of spatial nodes
Nt = int((tmax/ht + 0.5)) # number of time steps
Nx2 = int(Nx/2)

Psi = np.zeros(Nx+1, dtype=np.complex128)
Psi2 = np.zeros(Nx+1)
V = np.zeros(Nx+1)
x = np.zeros(Nx+1)

for i in range(1,Nx+1):
    x[i] = (i-Nx2-1)*hx
    V[i] = Pot(x[i], a, V0)

Init(Psi, x, Nx, x0, sig, k0)