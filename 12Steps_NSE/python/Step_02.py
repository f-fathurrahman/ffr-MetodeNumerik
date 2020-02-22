import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.style.use("dark_background")

Nx = 41
dx = 2.0/(Nx-1)
Nt = 20
dt = 0.025

u = np.ones(Nx)
# u = 2 between 0.5 and 1 as per IC
u[int(0.5/dx):int(1/dx+1)] = 2.0

x = np.linspace(0.0, 2.0, Nx)

plt.clf()
plt.plot(x, u)
plt.ylim(0.9, 2.1)
plt.savefig("IMG_Step_02_IC.png", dpi=150)

un = np.ones(Nx)

for n in range(Nt):
    un = u.copy()
    for i in range(1,Nx):
        u[i] = un[i] - un[i]*dt/dx*(un[i] - un[i-1])
    plt.clf()
    plt.plot(x, u)
    plt.ylim(0.9, 2.1)
    filename = "IMG_Step_02_{:04d}.png".format(n)
    plt.savefig(filename, dpi=150)

