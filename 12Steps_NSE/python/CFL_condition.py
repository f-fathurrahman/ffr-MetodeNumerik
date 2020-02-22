import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.style.use("dark_background")

def linear_convection(Nx):

    print("Nx = ", Nx)

    dx = 2.0/(Nx-1)
    Nt = 25
    dt = 0.025
    c = 1.0

    u = np.ones(Nx)
    # u = 2 between 0.5 and 1 as per IC
    u[int(0.5/dx):int(1/dx+1)] = 2.0

    x = np.linspace(0.0, 2.0, Nx)

    un = np.ones(Nx)

    n = 0
    plt.clf()
    plt.plot(x, u)
    plt.ylim(0.9, 2.1)
    plt.title("Nx = " + str(Nx))
    filename = "IMG_CFL_Nx_{:d}_{:04d}.png".format(Nx,n)
    plt.savefig(filename, dpi=150)

    for n in range(Nt):
        un = u.copy()
        for i in range(1,Nx):
            u[i] = un[i] - c*dt/dx*(un[i] - un[i-1])
        plt.clf()
        plt.plot(x, u)
        plt.ylim(0.9, 2.1)
        plt.title("Nx = " + str(Nx))
        filename = "IMG_CFL_Nx_{:d}_{:04d}.png".format(Nx,n+1)
        plt.savefig(filename, dpi=150)

linear_convection(41)
linear_convection(61)
linear_convection(71)
linear_convection(85)