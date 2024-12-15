import numpy as np

def solve_poisson2d_dirichlet( \
    f_func, g_func, \
    bx0, bxf, by0, byf, \
    D, Nx, Ny, \
    TOL=1e-8, NiterMax=500):

    # Get boundary information
    x0 = D[0]; xf = D[1]
    y0 = D[2]; yf = D[3]

    Δx = (xf - x0)/Nx
    Δy = (yf - y0)/Ny

    x = np.linspace(x0, xf, Nx+1)
    y = np.linspace(y0, yf, Ny+1)

    u = np.zeros( (Nx+1, Ny+1) )

    for j in range(Ny+1):
        u[0,j] = bx0( y[j] )
        u[Nx,j] = bxf( y[j] )

    for i in range(Nx+1):
        u[i,0] = by0( x[i] )
        u[i,Ny] = byf( x[i] )


    # Initial values for other nodes
    sum_of_bv = np.sum(u[0,:]) + np.sum(u[Nx,:]) + \
                np.sum(u[:,0]) + np.sum(u[:,Ny])
    u[1:Nx,1:Ny] = sum_of_bv / (2*(Nx+Ny-2))

    f = np.zeros( (Nx+1, Ny+1) )
    g = np.zeros( (Nx+1, Ny+1) )
    for i in range(Nx+1):
        for j in range(Ny+1):
            f[i,j] = f_func(x[i], y[j])
            g[i,j] = g_func(x[i], y[j])

    rxy = 0.5 * Δx**2 * Δy**2 / (Δx**2 + Δy**2)
    ry = 0.5 * Δy**2 / (Δx**2 + Δy**2)
    rx = 0.5 * Δx**2 / (Δx**2 + Δy**2)
    for iterJacobi in range(NiterMax):
        err = 0.0
        u_old = np.copy(u)
        # loop only for internal nodes
        for j in range(1,Ny):
            for i in range(1,Nx):
                u[i,j] = rx*(u[i,j+1] + u[i,j-1]) + ry*(u[i+1,j] + u[i-1,j]) + \
                         rxy*( g[i,j]*u[i,j] - f[i,j])
        # calculate error
        err = np.mean( (u - u_old)**2 )
        print("iterJacobi = %8d, err = %18.10e" % (iterJacobi, err))
        if err < TOL:
            print("Convergence is achieved")
            break
    return u, x, y


def f_func(x,y):
    return 0.0

def g_func(x,y):
    return 0.0

def bx0(y):
    return 0.0

def bxf(y):
    return 0.0

def by0(x):
    return np.sin(2*x)

def byf(x):
    return 0.0

Nx = 50
Ny = 50
D = [0.0, np.pi, 0.0, np.pi]

u, x, y = solve_poisson2d_dirichlet( \
    f_func, g_func, \
    bx0, bxf, by0, byf, \
    D, Nx, Ny, \
    TOL=1e-8, NiterMax=500)


X, Y = np.meshgrid(x, y)

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

fig = plt.figure()
plt.clf()
ax = fig.subplots(subplot_kw={"projection": "3d"})
plt.ylabel("y label")
plt.xlabel("x label")

# Plot the surface.
#surf = ax.plot_surface(X, Y, u, cmap=cm.coolwarm,
#                       linewidth=0, antialiased=True)
surf = ax.plot_surface(X, Y, u.T, linewidth=1, cmap=cm.coolwarm, antialiased=True)

# Customize the z axis.
#ax.set_zlim(-1.01, 1.01)
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter("%.02f"))

# Add a color bar which maps values to colors.
#fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
plt.savefig("TEMP_IMG_poisson2d.png", dpi=150)
