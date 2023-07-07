import numpy as np
import matplotlib.pyplot as plt
import scipy

plt.style.use("dark_background")


α = 1.0

# In case source is needed
# Here we simply set it it zero
def f_source(x, y, t):
    return np.zeros(np.size(x))


def create_D2_1d(N):
    D2 = scipy.sparse.lil_matrix((N,N))
    # first row
    D2[0,0:2] = np.array([-2.0, 1.0])
    # last row
    D2[N-1,N-2:N] = np.array([1.0, -2.0])
    # Other rows
    for irow in range(1,N-1):
        D2[irow,irow-1:irow+2] = np.array([1.0, -2.0, 1.0])
    return D2

def create_laplacian2d(Nx, Ny, dx, dy):
    D2x = create_D2_1d(Nx)/dx**2
    D2y = create_D2_1d(Ny)/dy**2

    Ix = scipy.sparse.diags(np.ones(Nx))
    Iy = scipy.sparse.diags(np.ones(Ny))

    # Laplacian matrix
    D2xy = scipy.sparse.kron(D2x, Iy) + scipy.sparse.kron(Ix, D2y)
    return D2xy





Lx = 1.0
Ly = 1.0

# Hollow square, centered at (Lx/2,Ly/2)
a = 0.25
b = 0.25
cx = Lx/2
cy = Ly/2

NsegX = 48
NsegY = 48

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

# Setup grid data (rectangular)
Npoints = Nx*Ny
rgrid = np.zeros((Npoints,2))
ip = 0
ip2xy = np.zeros((Npoints,2), dtype=np.int32)
xy2ip = np.zeros((Nx,Ny), dtype=np.int32)
for i in range(Nx):
    for j in range(Ny):
        rgrid[ip,0] = xgrid[i]
        rgrid[ip,1] = ygrid[j]
        ip2xy[ip,0] = i
        ip2xy[ip,1] = j
        xy2ip[i,j] = ip
        ip += 1


# Search for boundary points
OUTSIDE_BC_POINTS = []
ip = 0
for i in range(Nx):
    for j in range(Ny):
        if (i == 0) or (i == Nx-1):
            OUTSIDE_BC_POINTS.append(ip)
        elif (j == 0) or (j == Ny-1):
            OUTSIDE_BC_POINTS.append(ip)
        ip += 1


SMALL = 1e-10
IDX_BOUNDARY = []
IDX_HOLLOW = []

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
    elif is_boundary:
        IDX_BOUNDARY.append(ip)

T_inside = 0.0

u = np.zeros((Nx,Ny))

# Initial condition
u[:,:] = 15.0

# Set zero temp in hollow (or np.nan ?)
for ip in IDX_HOLLOW:
    i, j = ip2xy[ip,:]
    u[i,j] = 0.0



#
# Build Laplacian matrix
#
D2xy = scipy.sparse.lil_matrix(create_laplacian2d(Nx, Ny, Δx, Δy))

# Build b vector (zeros in case of Laplace equation)
b = np.zeros( Nx*Ny )
f = np.zeros( Nx*Ny )

# Needed for source evaluation and plotting
X, Y = np.meshgrid(xgrid, ygrid)

Ixy = scipy.sparse.diags(np.ones(Nx*Ny))

NtimeSteps = 100
t = 0.0
dt = 0.01

for ik in range(1,NtimeSteps+1):

    t = t + dt # Next time step

    print("Begin t = ", t)

    # Build LHS for implicit method
    LHS = scipy.sparse.lil_matrix( Ixy - α*D2xy*dt )

    # b-vector
    uflat = u.flatten()
    f[:] = f_source(X, Y, t).T.flatten()
    b[:] = uflat[:] + f*dt

    #
    # Apply BC for next time step
    #
    # Loop for all points (Ny) in the y-direction
    for j in range(Ny):
        u[0,j] = 100.0
        u[Nx-1,j] = 20.0
    # Loop for all points (Nx) in the x-direction
    for i in range(Nx):
        u[i,0] = 50.0
        u[i,Ny-1] = 20.0

    # Set u?

    # Apply BC to the Laplacian matrix and b vector
    uflat = u.flatten()
    for ip in OUTSIDE_BC_POINTS:
        LHS[ip,:] = np.array([0.0]*Nx*Ny)
        LHS[ip,ip] = 1.0
        b[ip] = uflat[ip]

    for ip in IDX_BOUNDARY:
        LHS[ip,:] = np.array([0.0]*Nx*Ny)
        LHS[ip,ip] = 1.0
        b[ip] = T_inside

    for ip in IDX_HOLLOW:
        LHS[ip,:] = np.array([0.0]*Nx*Ny)
        LHS[ip,ip] = 1.0
        b[ip] = np.nan # 0

    LHS = scipy.sparse.csr_matrix(LHS)
    u_next = scipy.sparse.linalg.spsolve(LHS, b)
    u_next = u_next.reshape(Nx, Ny)

    #du = u_next.flatten() - u.flatten()
    #print("norm du = ", np.dot(du.flatten(), du.flatten()))

    u[:,:] = u_next[:,:]

    print("t = ", t, " is done")


# Plot
"""
fig = plt.figure()
ax1 = fig.add_subplot(111)
str_title = "Solusi numerik t = {:10.5f}".format(t)
cnt = ax1.contourf(X, Y, u.T, cmap="coolwarm")
ax1.set_title(str_title)
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_aspect(1.0)
fig.colorbar(mappable=cnt)
plt.show()
"""

# Simpler way
plt.clf()
str_title = "Solusi numerik t = {:10.5f}".format(t)
plt.contourf(X, Y, u.T, cmap="coolwarm")
plt.title(str_title)
plt.xlabel("x")
plt.ylabel("y")
plt.gca().set_aspect(1.0)
plt.colorbar()
plt.show()
