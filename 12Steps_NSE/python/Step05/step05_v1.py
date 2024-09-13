from mpl_toolkits.mplot3d import Axes3D

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

Nx = 81
Ny = 81
Nt = 100
c = 1.0
Δx = 2.0/(Nx-1)
Δy = 2.0/(Ny-1)
σ = 0.2
Δt = σ*Δx

print("CFL = ", c*Δt/Δx)

x = np.linspace(0.0, 2.0, Nx)
y = np.linspace(0.0, 2.0, Ny)

u = np.ones( (Nx, Ny) )
un = np.ones( (Nx, Ny) )

# Initial conditions
X, Y = np.meshgrid(x,y)
idx_X = (X > 1.25) & (X < 1.75)
idx_Y = (Y > 0.5) & (Y < 1.0)
u[idx_X & idx_Y] = 2.0

plt.clf()
fig = plt.figure(figsize=(11,7))
ax = fig.add_subplot(projection="3d")
ax.plot_surface(X, Y, u[:], cmap=cm.viridis)
ax.set_zlim(1.0, 2.1)
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("IMG_Step_05_0000.png", dpi=150)


"""
for n in range(Nt + 1):
    #for n in range(2):
    un = u.copy()
    row, col = u.shape
    for j in range(1, row):
        for i in range(1, col):
            u[j, i] = (un[j, i] - (c*Δt/Δx * (un[j, i] - un[j, i - 1])) - (c*Δt/Δy * (un[j, i] - un[j - 1, i])))
            # fix the boundary conditions
            u[0, :] = 1.0
            u[-1, :] = 1.0
            u[:, 0] = 1.0
            u[:, -1] = 1.0
    ax.cla()
    ax.plot_surface(X, Y, u[:], cmap=cm.viridis)
    filename = "IMG_Step_05_{:04d}.png".format(n + 1)
    fig.savefig(filename, dpi=150)
    print("Done: ", n)
"""
