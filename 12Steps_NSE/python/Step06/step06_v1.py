from mpl_toolkits.mplot3d import Axes3D

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

Nx = 101
Ny = 101
Nt = 81
c = 1.0 # not used?
Δx = 2.0/(Nx-1)
Δy = 2.0/(Ny-1)
σ = 0.2
Δt = σ*Δx

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

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(projection="3d")
X, Y = np.meshgrid(x, y)

ax.plot_surface(X, Y, u, cmap=cm.viridis)
ax.set_xlabel("$x$")
ax.set_xlabel("$y$")
ax.set_zlim(1.0, 2.1)
plt.savefig("IMG_Step_06_0000.png", dpi=150)

from numba import jit

@jit(nopython=True)
def do_time_step(Δt, Δx, Δy, u, un, v, vn):
    Nx = u.shape[0]
    Ny = u.shape[1]
    for i in range(1,Nx):
        for j in range(1,Ny):
            u[i,j] = un[i,j] - \
                     un[i,j] * Δt/Δx * ( un[i,j] - un[i-1,j] ) - \
                     vn[i,j] * Δt/Δy * ( un[i,j] - un[i,j-1] )
            v[i,j] = vn[i,j] - \
                     un[i,j] * Δt/Δx * ( vn[i,j] - vn[i-1,j] ) - \
                     vn[i,j] * Δt/Δy * ( vn[i,j] - vn[i,j-1] )
    return


# warm-up
def do_warm_up():
    print("Doing warm up")
    Nx = 2
    Ny = 2
    u = np.zeros((Nx,Ny))
    v = np.zeros((Nx,Ny))
    un = np.copy(u)
    vn = np.copy(v)
    do_time_step(1.0, 1.0, 1.0, u, un, v, vn)
    print("Done doing warm up")
    return

do_warm_up()

import time
time_start = time.perf_counter()

for n in range(Nt + 1):
    #
    un = u.copy()
    vn = v.copy()
    #
    do_time_step(Δt, Δx, Δy, u, un, v, vn)
    #
    # Set BC
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

time_end = time.perf_counter()
print("Elapsed time: {} s".format(time_end - time_start))

ax.cla()
surf = ax.plot_surface(X, Y, u, cmap=cm.viridis)
ax.set_xlabel("$x$")
ax.set_xlabel("$y$")
ax.set_zlim(1.0, 2.1)
filename = "IMG_Step_06_last.png"
fig.savefig(filename, dpi=150)

