import numpy as np
import matplotlib.pyplot as plt

import scipy.sparse

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


def u_exact(x, y):
    return 1/np.sinh(2*np.pi) * np.sin(2*x) * np.sinh(2*(np.pi-y))

# Definisi syarat batas:
bx0 = lambda y: 0.0
bxf = lambda y: 0.0
by0 = lambda x: np.sin(2*x)
byf = lambda x: 0.0

x0 = 0.0; xf = np.pi
y0 = 0.0; yf = np.pi


Nx = 75; Ny = 75

dx = (xf-x0)/Nx
x = np.linspace(x0, xf, Nx+1)

dy = (yf - y0)/Ny
y = np.linspace(y0, yf, Ny+1)

Nx1 = Nx + 1
Ny1 = Ny + 1

u = np.zeros( (Nx+1,Ny+1) )

for j in range(Ny1):
    u[0,j] = bx0(y[j])
    u[Nx1-1,j] = bxf(y[j])

for i in range(Nx1):
    u[i,0] = by0(x[i])
    u[i,Ny-1] = byf(x[i])


D2x = create_D2_1d(Nx1)/dx**2
D2y = create_D2_1d(Ny1)/dy**2

Ix = scipy.sparse.diags(np.ones(Nx1))
Iy = scipy.sparse.diags(np.ones(Ny1))

# Laplacian matrix
# Still use lil_matrix because we need slicing operation (to apply BC)
D2xy = scipy.sparse.lil_matrix( scipy.sparse.kron(D2x, Iy) + scipy.sparse.kron(Ix, D2y) )
#D2xy = scipy.sparse.csr_matrix( scipy.sparse.kron(D2x, Iy) + scipy.sparse.kron(Ix, D2y) )


# Search for boundary points
bc_pts_list = []
ip = 0
for i in range(Nx1):
    for j in range(Ny1):
        if (i == 0) or (i == Nx):
            bc_pts_list.append(ip)
        elif (j == 0) or (j == Ny):
            bc_pts_list.append(ip)
        ip += 1

# Build b vector (zeros in case of Laplace equation)
b = np.zeros(Nx1*Ny1)

# Apply BC to the Laplacian matrix and b vector
uflat = u.flatten()
for ip in bc_pts_list:
    D2xy[ip,:] = np.array([0.0]*Nx1*Ny1)
    D2xy[ip,ip] = 1.0
    b[ip] = uflat[ip]

# Convert D2xy to CSR
D2xy = scipy.sparse.csr_matrix(D2xy)

u_sol = scipy.sparse.linalg.spsolve(D2xy, b)
u_sol_xy = u_sol.reshape(Nx1,Ny1) # convert to 2d array

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
X, Y = np.meshgrid(x, y)

du = u_exact(X, Y) - u_sol_xy.T
du_norm = np.dot(du.flatten(), du.flatten())/len(du)
print("du_norm = ", du_norm)

# Notice the transpose for u
ax.plot_surface(X, Y, u_sol_xy.T, cmap="coolwarm")

ax.set_title("Solusi numerik")
ax.set_xlabel("x")
ax.set_ylabel("y");
plt.show()
