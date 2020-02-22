import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.style.use("dark_background")

Nx = 41
dx = 2.0/(Nx-1)
Nt = 25
dt = 0.025
c = 1.0

u = np.ones(Nx)
# u = 2 between 0.5 and 1 as per IC
u[int(0.5/dx):int(1/dx+1)] = 2.0

x = np.linspace(0.0, 2.0, Nx)
plt.clf()
plt.plot(x, u)
plt.savefig("IMG_Step_01_IC.png", dpi=150)
