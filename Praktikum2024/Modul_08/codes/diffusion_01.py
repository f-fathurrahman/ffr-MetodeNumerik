# Explicit method (forward Euler) for diffusion equation
# Store solutions for all time steps

import numpy as np
import matplotlib.pyplot as plt

# Global variables !!!
L = 1.5
Tfinal = 1.0
α = 0.1

DO_PLOT = False

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

# Here we store the data for all time steps
u = np.zeros((Nt+1, Nx+1))

u[0,:] = initial_cond(x)
for n in range(0,Nt):
    # Apply boundary condition
    u[n,0] = 0.0 # syarat batas pada x=0
    u[n,Nx] = 0.0 # syarat batas pada x=L
    for i in range(1,Nx):
        fni = source_term(t[n], x[i])
        u[n+1,i] = u[n,i] + F*( u[n,i+1] - 2*u[n,i] + u[n,i-1] ) + Δt*fni

    if DO_PLOT:
        plt.clf()
        plt.plot(x, u[n,:])
        plt.title("t = " + str(t[n]))
        plt.ylim(0.0, 0.3)
        plt.savefig("IMG_diff1d_explicit_" + str(n) + ".png", dpi=150)
        print("n = " + str(n) + " is done")

if DO_PLOT:
    plt.clf()
    plt.plot(x, u[Nt,:], label="numerical")
    plt.plot(x, u_exact(t[Nt],x), label="exact sol")
    plt.title("t = " + str(t[Nt]))
    plt.ylim(0.0, 0.3)
    plt.legend()
    plt.savefig("IMG_diff1d_explicit_COMPARE_" + str(n) + ".png", dpi=150)

# Difference between exact solution and numerical solution
Δu = u_exact(t[Nt],x) - u[Nt,:]
norm_Δu = np.linalg.norm(Δu)
print("norm_du = ", norm_Δu)