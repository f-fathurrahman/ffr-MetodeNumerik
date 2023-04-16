import numpy as np
import matplotlib.pyplot as plt

def initial_cond(x, A=1.0, L=1.0, σ=0.02):
    return A*np.exp( -0.5*( (x - L/10)/σ )**2 )


def initial_cond_v2(x, A=1.0, L=1.0, σ=0.02):
    if x < L/5:
        return A*cos(5*pi/L*(x-L/10))
    else:
        return 0.0

L = 1.0
v = 1.0
A = 1.0
σ = 0.02

Δx = 0.01
x = np.arange(0.0, 1.0, Δx)

# Initial condition
u = initial_cond(x, A=A, L=L, σ=σ)
#u = initial_cond_v2.(x; A=A, L=L, σ=σ)

Δt = 0.001
C = v*Δt/Δx

Nx = len(x)
unp1 = np.zeros(Nx)

plt.clf()
plt.plot(x, u, label="u")
plt.ylim(-0.1, 1.1)
plt.grid(True)
plt.legend()
plt.savefig("IMG_adv_" + str(0) + ".png", dpi=150)

INTERVAL_DO_PLOT = 50

for n in range(1,1001): # loop over time
    print("Time = %.5f" % (n*Δt))
    unp1[0] = 0.0 # Boundary condition
    unp1[Nx-1] = 0.0 # Boundary condition
    for i in range(1,Nx-1):
        unp1[i] = u[i] - 0.5*C*( u[i+1] - u[i-1] )
    u[:] = unp1[:] # copy new value to old value
    
    if (n % INTERVAL_DO_PLOT == 0):
        print("Plotting the function")
        plt.clf()
        plt.plot(x, u, label="u")
        plt.ylim(-0.1, 1.1)
        plt.grid(True)
        plt.legend()
        plt.savefig("IMG_adv_" + str(n) + ".png", dpi=150)

