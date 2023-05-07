import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse

Nx = 3
Ny = 4

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

D2x = create_D2_1d(Nx)
D2y = create_D2_1d(Ny)

Ix = scipy.sparse.diags(np.ones(Nx))
Iy = scipy.sparse.diags(np.ones(Ny))

D2xy = scipy.sparse.lil_matrix( scipy.sparse.kron(D2x, Iy) + scipy.sparse.kron(Ix, D2y) )

plt.imshow(D2xy.todense(), interpolation="none")
plt.show()

bc_pts_list = []
ip = 0
for i in range(Nx):
    for j in range(Ny):
        if (i == 0) or (i == Nx-1):
            bc_pts_list.append(ip)
        elif (j == 0) or (j == Ny-1):
            bc_pts_list.append(ip)
        ip += 1

for ip in bc_pts_list:
    D2xy[ip,:] = np.array([0.0]*Nx*Ny)
    D2xy[ip,ip] = 1.0

