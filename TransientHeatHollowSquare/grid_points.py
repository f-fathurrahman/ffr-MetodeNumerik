import numpy as np
import matplotlib.pyplot as plt

plt.style.use("dark_background")

Lx = 1.0
Ly = 1.0

# Hollow square, centered at (Lx/2,Ly/2)
a = 0.5
b = 0.5
cx = Lx/2
cy = Ly/2

NsegX = 8
NsegY = 8

# Current limitation
assert NsegX % 4 == 0
assert NsegY % 4 == 0

Nx = NsegX + 1
Ny = NsegY + 1

xgrid = np.linspace(0.0, Lx, Nx)
Δx = xgrid[1] - xgrid[0]

ygrid = np.linspace(0.0, Ly, Ny)
Δy = ygrid[1] - ygrid[0]

print("Δx = ", Δx)
print("Δy = ", Δy)

print("xgrid = ", xgrid)

Npoints = Nx*Ny
rgrid = np.zeros((Npoints,2))
ip = 0
for i in range(Nx):
    for j in range(Ny):
        rgrid[ip,0] = xgrid[i]
        rgrid[ip,1] = ygrid[j]
        ip += 1

SMALL = 1e-10

IDX_BOUNDARY = []
IDX_HOLLOW = []

plt.clf()
for ip in range(Npoints):
    x = rgrid[ip,0]
    y = rgrid[ip,1]
    
    is_x = abs(x - cx) < a/2
    is_y = abs(y - cy) < b/2
    in_hollow = is_x and is_y

    is_boundary_x = (abs(x - cx) - a/2) <= SMALL
    is_boundary_y = (abs(y - cy) - b/2) <= SMALL
    is_boundary = is_boundary_x and is_boundary_y

    if in_hollow:
        IDX_HOLLOW.append(ip) 
        plt.plot([x], [y], marker="o", color="gray")
    elif is_boundary:
        IDX_BOUNDARY.append(ip)
        plt.plot([x], [y], marker="o", color="green")
    else:
        plt.plot([x], [y], marker="o", color="red")


#for p in internal_corner_points:
#    plt.plot([p[0]], [p[1]], marker="*", color="magenta")

internal_corner_points = [
    [0.5 - a/2, 0.5 - b/2],
    [0.5 + a/2, 0.5 - b/2],
    [0.5 - a/2, 0.5 + b/2],
    [0.5 + a/2, 0.5 + b/2]
]

plt.axis("equal")
plt.savefig("IMG_grid_points.png", dpi=150)

