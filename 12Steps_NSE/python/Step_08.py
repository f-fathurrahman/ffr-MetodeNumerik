from mpl_toolkits.mplot3d import Axes3D

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import matplotlib
matplotlib.style.use("dark_background")

Nx = 41
Ny = 41
Nt = 120
c = 1.0
Δx = 2.0/(Nx-1)
Δy = 2.0/(Ny-1)
σ = 0.0009
ν = 0.01
Δt = σ*Δx*Δy/ν

x = np.linspace(0.0, 2.0, Nx)
y = np.linspace(0.0, 2.0, Ny)

u = np.ones( (Nx,Ny) )
v = np.ones( (Nx,Ny) )

un = np.ones( (Nx,Ny) )
vn = np.ones( (Nx,Ny) )

# Initial condition
# set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
idx_xi = int(0.5/Δx)
idx_xf = int(1.0/Δx + 1)
idx_yi = int(0.5/Δy)
idx_yf = int(1.0/Δy + 1)

u[idx_xi:idx_xf, idx_yi:idx_yf] = 2.0
v[idx_xi:idx_xf, idx_yi:idx_yf] = 2.0

fig = plt.figure(figsize=(11,7))
ax = fig.gca(projection="3d")
X, Y = np.meshgrid(x, y)

ax.plot_surface(X, Y, u, cmap=cm.jet)
ax.set_xlabel("$x$")
ax.set_xlabel("$y$")
ax.set_zlim(1.0, 2.1)
plt.savefig("IMG_Step_08_0000.png", dpi=150)

for n in range(Nt + 1):
    #
    un = u.copy()
    vn = v.copy()
    #
    u[1:-1, 1:-1] = ( un[1:-1, 1:-1] - 
                      Δt/Δx * un[1:-1, 1:-1] * ( un[1:-1, 1:-1] - un[1:-1, 0:-2] ) -
                      Δt/Δy * vn[1:-1, 1:-1] * ( un[1:-1, 1:-1] - un[0:-2, 1:-1] ) +
                      ν*Δt/Δx**2 * ( un[1:-1,2:] - 2*un[1:-1, 1:-1] + un[1:-1, 0:-2] ) + 
                      ν*Δt/Δy**2 * ( un[2:, 1:-1] - 2*un[1:-1, 1:-1] + un[0:-2, 1:-1] )
                    )

    v[1:-1, 1:-1] = ( vn[1:-1, 1:-1] -
                      Δt/Δx * un[1:-1, 1:-1] * ( vn[1:-1, 1:-1] - vn[1:-1, 0:-2] ) -
                      Δt/Δy * vn[1:-1, 1:-1] * ( vn[1:-1, 1:-1] - vn[0:-2, 1:-1] ) +
                      ν*Δt/Δx**2 * ( vn[1:-1, 2:] - 2*vn[1:-1, 1:-1] + vn[1:-1, 0:-2] ) +
                      ν*Δt/Δy**2 * ( vn[2:, 1:-1] - 2*vn[1:-1, 1:-1] + vn[0:-2, 1:-1] )
                    )
    #
    u[0,  :] = 1.0
    u[-1, :] = 1.0
    u[:,  0] = 1.0
    u[:, -1] = 1.0
    #
    v[0,  :] = 1.0
    v[-1, :] = 1.0
    v[:,  0] = 1.0
    v[:, -1] = 1.0

    fig.clf()
    ax = fig.gca(projection="3d")
    surf = ax.plot_surface(X, Y, u, cmap=cm.jet)
    ax.set_xlabel("$x$")
    ax.set_xlabel("$y$")
    ax.set_zlim(1.0, 2.1)
    filename = "IMG_Step_08_{:04d}.png".format(n + 1)
    fig.savefig(filename, dpi=150)
    print("Done: ", n)