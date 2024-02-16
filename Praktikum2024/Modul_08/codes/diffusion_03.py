# Implicit method (forward Euler) for diffusion equation
# Using two vectors for current and future time

import numpy as np
import matplotlib.pyplot as plt

# Global variables !!!
L = 1.5
Tfinal = 1.0
α = 1.0

DO_PLOT = True

# Manufactured solution
def u_exact(t, x):
    return 5*x*t*(L-x)

def initial_cond(x):
    return u_exact(0, x)

def source_term(t, x):
    return 10*α*t + 5*x*(L - x)

Nx = 25
x = np.linspace(0.0, L, Nx+1)
Δx = x[1] - x[0]
print("Δx = ", Δx)

Nt = 200
t = np.linspace(0.0, Tfinal, Nt+1)
Δt = t[1] - t[0]
print("Δt = ", Δt)
print("Final t = ", t[-1])

F = α * Δt / Δx**2
print("F = ", F)


# Build matrix A
A = np.zeros( (Nx-1,Nx-1) )

A[0,0] = 1 + 2*F
A[0,1] = -F
for i in range(1,Nx-2):
    A[i,i] = 1 + 2*F
    A[i,i+1] = -F
    A[i,i-1] = -F
A[Nx-2,Nx-2] = 1 + 2*F
A[Nx-2,Nx-3] = -F

# Here we store the data for all time steps
un = np.zeros(Nx+1)
unp1 = np.zeros(Nx+1)

# time loop
b = np.zeros(Nx-1)
for n in range(1,Nt+1):
    for i in range(Nx-1):
        b[i] = un[i+1] + Δt*source_term(t[n], x[i+1])
    # Selesaikan sistem persamaan linear
    unp1[1:Nx] = np.linalg.solve(A, b)
    #
    un[:] = unp1[:] # for next iteration

if DO_PLOT:
    plt.clf()
    plt.plot(x, un, label="numerical")
    plt.plot(x, u_exact(t[Nt],x), label="exact sol")
    plt.title("t = " + str(t[Nt]))
    plt.legend()
    plt.savefig("IMG_diffusion1d_implicit_COMPARE_" + str(n) + ".png", dpi=150)

# Difference between exact solution and numerical solution
Δu = u_exact(t[Nt],x) - un
norm_Δu = np.linalg.norm(Δu)
print("norm_du = ", norm_Δu)
