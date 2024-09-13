from mpl_toolkits.mplot3d import Axes3D

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

Nx = 31
Ny = 31
Nt = 50
ν = 0.05
Δx = 2.0/(Nx-1)
Δy = 2.0/(Ny-1)
σ = 0.25
Δt = σ*Δx*Δy/ν

x = np.linspace(0.0, 2.0, Nx)
y = np.linspace(0.0, 2.0, Ny)

u = np.ones( (Nx,Ny) )
un = np.ones( (Nx,Ny) )

# Initial condition
# set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
iΔx_xi = int(0.5/Δx)
iΔx_xf = int(1.0/Δx + 1)
iΔx_yi = int(0.5/Δy)
iΔx_yf = int(1.0/Δy + 1)

u[iΔx_xi:iΔx_xf, iΔx_yi:iΔx_yf] = 2.0

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(projection="3d")
X, Y = np.meshgrid(x, y)

ax.plot_surface(X, Y, u, cmap=cm.viridis)
ax.set_xlabel("$x$")
ax.set_xlabel("$y$")
ax.set_xlim(0, 2)
ax.set_ylim(0, 2)
ax.set_zlim(1, 2.5)

plt.savefig("IMG_Step_07_0000.png", dpi=150)


from numba import jit

@jit(nopython=True)
def do_time_step(u, un, ν, Δt, Δx, Δy):
    Nx = u.shape[0]
    Ny = u.shape[1]
    for i in range(1,Nx-1):
        for j in range(1,Ny-1):
            u[i,j] = un[i,j] + \
                     ν*Δt/Δx**2 * ( un[i+1,j] - 2*un[i,j] + un[i-1,j] ) + \
                     ν*Δt/Δy**2 * ( un[i,j+1] - 2*un[i,j] + un[i,j-1] )
    return


# warm-up
def do_warm_up():
    print("Doing warm up")
    Nx = 2
    Ny = 2
    u = np.zeros((Nx,Ny))
    un = np.copy(u)
    do_time_step(u, un, 1.0, 1.0, 1.0, 1.0)
    print("Done doing warm up")
    return


do_warm_up()


import time
time_start = time.perf_counter()

for n in range(Nt + 1):
    un = u.copy()
    #
    do_time_step(u, un, ν, Δt, Δx, Δy)
    # set BC
    u[0,  :] = 1
    u[-1, :] = 1
    u[:,  0] = 1
    u[:, -1] = 1

time_end = time.perf_counter()
print("Elapsed time: {} s".format(time_end - time_start))

ax.cla()
ax.plot_surface(X, Y, u, cmap=cm.viridis)
ax.set_xlabel("$x$")
ax.set_xlabel("$y$")
ax.set_xlim(0, 2)
ax.set_ylim(0, 2)
ax.set_zlim(1, 2.5)
filename = "IMG_Step_07_last.png"
fig.savefig(filename, dpi=150)
