import numpy as np
import scipy.sparse
import matplotlib.pyplot as plt

α = 1.0

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

Δx = 1.0
Δy = 0.5
Nx = 3
Ny = 4

Laplacian = create_laplacian2d(Nx, Ny, Δx, Δy)

Lx = (Nx-1)*Δx
Ly = (Ny-1)*Δy

# Setup grid data (rectangular)
xgrid = np.linspace(0.0, Lx, Nx)
Δx = xgrid[1] - xgrid[0]

ygrid = np.linspace(0.0, Ly, Ny)
Δy = ygrid[1] - ygrid[0]

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


A = scipy.sparse.lil_matrix((Npoints,Npoints))
denom = Δx**2 * Δy**2
for ip in range(Npoints):
    
    print("")
    print("Row ip = ", ip)
    A[ip,ip] += -2*(Δx**2 + Δy**2)/denom

    i, j = ip2xy[ip,:]
    ii = i + 1
    if (ii >= 0) and (ii < Nx):
        ip1 = xy2ip[ii,j]
        print("col_idx = ", ip1)
        A[ip,ip1] += Δy**2/denom

    ii = i - 1
    if (ii >= 0) and (ii < Nx):
        ip1 = xy2ip[ii,j]
        print("col_idx = ", ip1)
        A[ip,ip1] += Δy**2/denom

    jj = j + 1
    if (jj >= 0) and (jj < Ny):
        ip1 = xy2ip[i,jj]
        print("col_idx = ", ip1)
        A[ip,ip1] += Δx**2/denom

    jj = j - 1
    if (jj >= 0) and (jj < Ny):
        ip1 = xy2ip[i,jj]
        print("col_idx = ", ip1)
        A[ip,ip1] += Δx**2/denom

