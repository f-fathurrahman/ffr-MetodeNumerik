import numpy as np
from math import *

def poisson2d_dirichlet(u, x, y, Nx, Ny, TOL, func):
    """
    Nx, Ny: no of nodes in x and y direction
    """
    MaxIter = 10000

    f = np.zeros((Nx,Ny))

    hx = (x[Nx-1] - x[0])/(Nx-1)
    kx = 1.0/(hx*hx)

    hy = (x[Ny-1] - y[0])/(Nx-1)
    ky = 1.0/(hy*hy)

    kxy = 2.0*(kx + ky)

    # calculate array for RHS
    for j in range(1,Ny-1):
        for i in range(1,Nx-1):
            f[i,j] = func(x[i],y[j])

    for iter in range(1,MaxIter+1):
        err = 0.0
        u_old = np.copy(u)
        # loop only for internal nodes
        for j in range(1,Ny-1):
            for i in range(1,Nx-1):
                u[i,j] = ( kx*( u[i-1,j] + u[i+1,j] ) +
                        ky*( u[i,j-1] + u[i,j+1] ) - f[i,j] ) / kxy
        # calculate error
        err = np.sum(np.abs(u - u_old))
        print('iter = %8d, err = %18.10e' % (iter,err))
        if err < TOL:
            print('Convergence is achieved')
            break


    return u


def func(x,y):
    return sin(x-y) - cos(x-y)

Nx = 21
Ny = 21
x = np.linspace(-pi,pi,Nx)
y = np.linspace(-pi,pi,Ny)

u0 = np.zeros((Nx,Ny))
# set BC
u0[0,:] = 0.0
u0[:,0] = 0.0
u0[Nx-1,:] = 0.0
u0[:,Ny-1] = 0.0

u = poisson2d_dirichlet(u0,x,y,Nx,Ny,1e-5,func)

X, Y = np.meshgrid(x, y)

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

fig = plt.figure()
plt.clf()
ax = fig.gca(projection='3d')


# Plot the surface.
#surf = ax.plot_surface(X, Y, u, cmap=cm.coolwarm,
#                       linewidth=0, antialiased=True)
surf = ax.plot_surface(X, Y, u, linewidth=1, cmap=cm.jet, antialiased=True)

# Customize the z axis.
#ax.set_zlim(-1.01, 1.01)
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
#fig.colorbar(surf, shrink=0.5, aspect=5)

plt.savefig('poisson2d_v4.png',dpi=200)

plt.show()

