import numpy as np
import matplotlib.pyplot as plt

L = 1.0
T = 0.05
alpha = 1.0
sigma = 0.05

def initial_cond(x):
    return np.exp( -(x-0.5*L)**2 / (2*sigma**2) )

def source_term(t, x):
    return 0.0

Nx = 50
x = np.linspace(0.0, L, Nx+1)
dx = x[1] - x[0]
print("dx = ", dx)

Nt = 300
t = np.linspace(0.0, T, Nt+1)
dt = t[1] - t[0]
print("dt = ", dt)

F = alpha * dt / dx**2
print("F = ", F)
if F > 0.5:
    print("WARNING: solution is not stable")

u = np.zeros((Nt+1, Nx+1))
for i in range(0,Nx+1):
    u[0,i] = initial_cond(x[i])

for n in range(0,Nt):
    u[n,0] = 0.0 # syarat batas pada x=0
    u[n,Nx] = 0.0 # syarat batas pada x=L
    for i in range(1,Nx):
        fni = source_term(t[n], x[i])
        u[n+1,i] = u[n,i] + F*( u[n,i+1] - 2*u[n,i] + u[n,i-1] ) + dt *  fni

    plt.clf()
    plt.plot(x, u[n,:])
    plt.title("t = " + str(t[n]))
    plt.ylim(0.0, 1.1)
    plt.savefig("IMG_diff1d_" + str(n) + ".png", dpi=150)

    print("n = " + str(n) + " is done")
