import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.style.use("dark_background")

Nx = 41
dx = 2.0/(Nx-1)
Nt = 20
nu = 0.3
sigma = 0.2
dt = sigma*dx**2/nu

x = np.linspace(0.0, 2.0, Nx)
u = np.ones(Nx)
u[int(0.5/dx):int(1.0/dx+1)] = 2.0

un = np.ones(Nx)

n = 0
plt.clf()
plt.plot(x, u)
plt.ylim(0.9, 2.1)
filename = "IMG_Step_03_{:04d}.png".format(0)
plt.savefig(filename, dpi=150)

for n in range(Nt):
    un = u.copy()
    for i in range(1,Nx-1):
        u[i] = un[i] + nu*dt/dx**2 * (un[i+1] - 2*un[i] + un[i-1])
    plt.clf()
    plt.plot(x, u)
    plt.ylim(0.9, 2.1)
    filename = "IMG_Step_03_{:04d}.png".format(n+1)
    plt.savefig(filename, dpi=150)

