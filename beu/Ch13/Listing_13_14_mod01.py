# Reflexion/transmission of quantum wave packet using tridiagonal solver
from math import pi, sqrt, fabs
import numpy as np
import cmath

import matplotlib.pyplot as plt

import matplotlib
matplotlib.style.use("dark_background")

import matplotlib_inline
matplotlib_inline.backend_inline.set_matplotlib_formats("svg")


def linsolve_tridiag(N, a, b, c, d):
    # Solves a system with tridiagonal matrix by LU factorization (diag(L) = 1).
    # a - lower codiagonal (i=2,n)
    # b - main diagonal (i=1,n)
    # c - upper codiagonal (i=1,n-1)
    # d - constant terms (i=1,n); solution on exit
    # n - order of system.
    if b[1] == 0.0:
        raise RuntimeError("TriDiagSys: singular matrix !")

    for i in range(2,N+1):
        # factorization
        a[i] = a[i]/b[i-1]
        b[i] = b[i] - a[i]*c[i-1]
        if b[i] == 0.0:
            raise RuntimeError("TriDiagSys: singular matrix !")
        d[i] = d[i] - a[i]*d[i-1]
    
    d[N] = d[N]/b[N]
    # backward substitution
    for i in range(N-1,0,-1):
        d[i] = (d[i] - c[i]*d[i+1])/b[i]
    return

def propagate(Psi, V, Nx, Δx, Δt):
    # Propagates the solution PSI = Psi + i Chi of the 1D Schrodinger equation
    #
    # i d/dt PSI(x,t) = [-(1/2) d2/dx2 - V(x)] PSI(x,t),
    # over the time interval Δt.
    # Uses the Crank-Nicolson scheme on a grid with
    # nx nodes and spacing hx and solves the tridiagonal discretized system by
    # LU factorization. Uses complex arithmetic.
    a = np.zeros(Nx+1, dtype=np.complex128) # diagonals of discretized system
    b = np.zeros(Nx+1, dtype=np.complex128)
    c = np.zeros(Nx+1, dtype=np.complex128)
    d = np.zeros(Nx+1, dtype=np.complex128)

    λ = Δt/(4.0*Δx*Δx)
    b[1] = 1.0
    c[1] = 0.0
    # build coefficients of discretized system
    for i in range(2,Nx):
        W = 2*λ + 0.5*Δt*V[i]
        a[i] = -λ
        b[i] = W - 1j
        c[i] = -λ
        d[i] = λ*Psi[i-1] - (W + 1j)*Psi[i] + λ*Psi[i+1] # constant term
    
    a[Nx] = 0.0
    b[Nx] = 1.0
    d[Nx] = 0.0
    # solve tridiagonal discretized system
    linsolve_tridiag(Nx, a, b, c, d)
    return d


def calc_potential_square(x, a, V0):
    if fabs(x) <= 0.5*a:
        return V0
    else:
        return 0

# Initial Gaussian wave packet
def initial_wavepacket(Psi, x, nx, x0, sig, k0):
    # Psi(x,0) = 1/sqrt(sqrt(2*pi)*sig) * exp[-(x-x0)^2/(4*sig^2)] * exp(ikx)
    # x0 - position of center, sig - half-width, k0 - average wave number
    a =  1.0/sqrt(sqrt(2*pi)*sig)
    b = -1.0/(4*sig*sig)
    for i in range(1,nx+1):
        dx = x[i] - x0
        f = a * np.exp(b*dx*dx)
        if f < 1e-10:
            f = 0.0
        Psi[i] = f*cmath.exp(1j*k0*x[i]) # probably np.math is also can be used


# Calculates the probability density Psi2[] of the wave function Psi[]
def calc_prob_dens(Psi, Psi2, nx, hx):
    for i in range(1,nx+1):
        Psi2[i] = abs(Psi[i])*abs(Psi[i])
        # unnormalized probability density
        if (Psi2[i] <= 1e-10):
            Psi2[i] = 0.0

    # integral by trapezoidal rule    
    PsiNorm = 0.5*(Psi2[1] + Psi2[nx])
    for i in range(2,nx):
        PsiNorm = PsiNorm + Psi2[i]
    PsiNorm = PsiNorm*hx

    for i in range(1,nx+1):
        Psi2[i] = Psi2[i]/PsiNorm

    return PsiNorm # normalized prob. density

# main
a = 5.0     # width of potential barrier
V0 = 30.0   # height of potential barrier
x0 = -20.0  # initial position of wave packet
sig = 1.0    # half-width of packet
k0 = 5.0   # average wave number of packet
xmax = 100.0
Δx = 5e-2 # spatial step size
tmax = 5.0
Δt = 5e-3 # time step
nout = 50 # output every nout steps

Nx = 2*int((xmax/Δx + 0.5) + 1) # odd number of spatial nodes
Nt = int((tmax/Δt + 0.5)) # number of time steps
Nx2 = int(Nx/2)
print("Nx = ", Nx)
print("Nt = ", Nt)

Psi = np.zeros(Nx+1, dtype=np.complex128)
Psi2 = np.zeros(Nx+1)
V = np.zeros(Nx+1)
x = np.zeros(Nx+1)
for i in range(1,Nx+1):
    x[i] = (i - Nx2 - 1)*Δx
    #V[i] = calc_potential_square(x[i], a, V0)

initial_wavepacket(Psi, x, Nx, x0, sig, k0)
Psir = np.real(Psi)
Psii = np.imag(Psi)
plt.figure(figsize=(6,3))
plt.plot(x, Psir, label="real")
plt.plot(x, Psii, label="imag")
plt.xlim(x0-5,x0+5)
plt.legend()
plt.savefig("IMG_initial.png", dpi=150);

for it in range(1,Nt+1):
    t = it*Δt
    print("time = %.5f" % t)
    Psi = propagate(Psi, V, Nx, Δx, Δt) # propagate solution by tridiagonal solver
    PsiNorm = calc_prob_dens(Psi, Psi2, Nx, Δx) # probability density
    if (it % nout == 0 or it == Nt):
        # output every nout steps
        fname = "IMG_psi2_{:08d}.png".format(it)
        plt.clf()
        plt.figure(figsize=(10,3))
        plt.plot(x[1:Nx+1], V[1:Nx+1]/10)
        plt.plot(x[1:Nx+1], Psi2[1:Nx+1])
        plt.xlim(-25.0, 25.0)
        plt.ylim(0.0, 0.5)
        plt.savefig(fname, dpi=150)
        plt.close()

