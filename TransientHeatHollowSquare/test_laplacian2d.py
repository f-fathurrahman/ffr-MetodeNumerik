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
