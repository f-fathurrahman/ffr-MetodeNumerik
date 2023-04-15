# Explicit method (forward Euler) for diffusion equation
# Using two vectors for current and future time

import numpy as np
import matplotlib.pyplot as plt

# Global variables !!!
L = 1.5
Tfinal = 1.0
α = 0.1

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
if F > 0.5:
    print("WARNING: solution is not stable")

# exit()

# Use only two vectors for the solution
un = np.zeros(Nx+1)
unp1 = np.zeros(Nx+1)

un[:] = initial_cond(x)
for n in range(0,Nt):
    # Apply boundary condition
    un[0] = 0.0 # syarat batas pada x=0
    un[Nx] = 0.0 # syarat batas pada x=L
    for i in range(1,Nx):
        fni = source_term(t[n], x[i])
        unp1[i] = un[i] + F*( un[i+1] - 2*un[i] + un[i-1] ) + Δt*fni
    #
    un[:]  = unp1[:] # update for the next iteration
    if DO_PLOT:
        plt.clf()
        plt.plot(x, un)
        plt.title("t = " + str(t[n]))
        plt.savefig("IMG_diff1d_explicit_" + str(n) + ".png", dpi=150)
        print("n = " + str(n) + " is done")

if DO_PLOT:
    plt.clf()
    plt.plot(x, un, label="numerical")
    plt.plot(x, u_exact(t[Nt],x), label="exact sol")
    plt.title("t = " + str(t[Nt]))
    plt.legend()
    plt.savefig("IMG_diff1d_explicit_COMPARE_" + str(n) + ".png", dpi=150)

# Difference between exact solution and numerical solution
Δu = u_exact(t[Nt],x) - un
norm_Δu = np.linalg.norm(Δu)
print("norm_du = ", norm_Δu)
