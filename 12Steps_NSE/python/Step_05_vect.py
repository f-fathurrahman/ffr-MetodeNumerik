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

x = np.linspace(0.0, 2.0, Nx)
y = np.linspace(0.0, 2.0, Ny)

u = np.ones( (Nx,Ny) )
un = np.ones( (Nx,Ny) )

# Initial conditions
u[int(0.5/Δy):int(1.0/Δy + 1),int(0.5/Δx):int(1.0/Δx + 1)] = 2.0

plt.clf()

fig = plt.figure(figsize=(11,7))
ax = fig.gca(projection="3d")
X, Y = np.meshgrid(x,y)

ax.plot_surface(X, Y, u[:], cmap=cm.jet)
ax.set_zlim(1.0, 2.1)

plt.savefig("IMG_Step_05_0000.png", dpi=150)

for n in range(Nt + 1):
    #for n in range(2):
    un = u.copy()
    u[1:, 1:] = (un[1:, 1:] - (c*Δt/Δx * (un[1:, 1:] - un[1:, :-1])) - (c*Δt/Δy * (un[1:, 1:] - un[:-1, 1:])))

    #fig.clf()
    #ax = fig.gca(projection="3d")
    #surf = ax.plot_surface(X, Y, u[:], cmap=cm.jet)
    #ax.set_zlim(1.0, 2.1)
    #filename = "IMG_Step_05_{:04d}.png".format(n + 1)
    #fig.savefig(filename, dpi=150)
    #print("Done: ", n)

