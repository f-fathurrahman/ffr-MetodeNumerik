import numpy as np
import scipy.sparse
import matplotlib.pyplot as plt

α = np.sqrt(1e-2)

# From Yang - Applied Numerical Methods with MATLAB 
# Example 9.3

def initial_cond(X,Y):
    return np.zeros(X.shape)

def bxyt(x,y,t):
    return np.exp(y)*np.cos(x) - np.exp(x)*np.cos(y)

bx0 = lambda y, t, x0: bxyt(x0, y, t)
bxf = lambda y, t, xf: bxyt(xf, y, t)
by0 = lambda x, t, y0: bxyt(x, y0, t)
byf = lambda x, y, yf: bxyt(x, yf, t)

def f_source(x, y, t):
    return np.zeros(x.shape)

# Build 2nd derivative matrix in 1d
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

# Build 2nd derivative matrix in 1d
def create_laplacian2d(Nx, Ny, dx, dy):
    D2x = create_D2_1d(Nx)/dx**2
    D2y = create_D2_1d(Ny)/dy**2

    Ix = scipy.sparse.diags(np.ones(Nx))
    Iy = scipy.sparse.diags(np.ones(Ny))

    # Laplacian matrix
    D2xy = scipy.sparse.kron(D2x, Iy) + scipy.sparse.kron(Ix, D2y)
    return D2xy


x0 = 0.0
xf = 4.0
y0 = 0.0
yf = 4.0

t0 = 0
tf = 5000
Nx = 40
Ny = 40
NtimeSteps = 50

Nx1 = Nx + 1
Ny1 = Ny + 1

#
dt = tf/NtimeSteps
dx = (xf - x0)/Nx
dy = (yf - y0)/Ny

x = np.linspace(x0, xf, Nx+1)
y = np.linspace(y0, yf, Ny+1)
t = np.linspace(t0, tf, NtimeSteps+1)

X, Y = np.meshgrid(x, y)


u = np.zeros((Nx+1,Ny+1))

# Search for boundary points
bc_pts_list = []
ip = 0
for i in range(Nx+1):
    for j in range(Ny+1):
        if (i == 0) or (i == Nx):
            bc_pts_list.append(ip)
        elif (j == 0) or (j == Ny):
            bc_pts_list.append(ip)
        ip += 1

# Initial condition
u[:,:] = initial_cond(X,Y).T

# Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
X, Y = np.meshgrid(x, y)
ax.plot_surface(X, Y, u.T, cmap="coolwarm")
ax.set_title("Initial")
ax.set_xlabel("x")
ax.set_ylabel("y")
plt.show()

D2xy = scipy.sparse.lil_matrix(create_laplacian2d(Nx1, Ny1, dx, dy))

# Build b vector (zeros in case of Laplace equation)
b = np.zeros( Nx1*Ny1 )
f = np.zeros( Nx1*Ny1 )

Ixy = scipy.sparse.diags(np.ones(Nx1*Ny1))

t = 0.0
for ik in range(1,NtimeSteps+1):

    t = t + dt # Next time step

    print("Begin t = ", t)

    # Build LHS
    LHS = scipy.sparse.lil_matrix( Ixy - α*D2xy*dt )

    # b-vector
    uflat = u.flatten()
    f[:] = f_source(X, Y, t).T.flatten()
    b[:] = uflat[:] + f*dt # previous u, flattened

    #
    # Apply BC for next time step
    #
    # Loop for all points (Ny+1) in the y-direction
    for j in range(Ny1):
        u[0,j] = bx0(y[j], t, x0)
        u[Nx,j] = bxf(y[j], t, xf)
    # Loop for all points (Nx+1) in the x-direction
    for i in range(Nx1):
        u[i,0] = by0(x[i], t, y0)
        u[i,Ny] = byf(x[i], t, yf)

    # Apply BC to the Laplacian matrix and b vector
    uflat = u.flatten()
    for ip in bc_pts_list:
        LHS[ip,:] = np.array([0.0]*Nx1*Ny1)
        LHS[ip,ip] = 1.0
        b[ip] = uflat[ip]

    LHS = scipy.sparse.csr_matrix(LHS)
    u_next = scipy.sparse.linalg.spsolve(LHS, b)
    u_next = u_next.reshape(Nx1, Ny1)

    du = u_next.flatten() - u.flatten()
    print("norm du = ", np.dot(du.flatten(), du.flatten()))

    u[:,:] = u_next[:,:]

    print("t = ", t, " is done")


# Plot
fig = plt.figure()

ax1 = fig.add_subplot(111, projection="3d")
str_title = "Solusi numerik t = {:10.5f}".format(t)
ax1.plot_surface(X, Y, u.T, cmap="coolwarm",  edgecolors="gray")
ax1.set_title(str_title)
ax1.set_xlabel("x")
ax1.set_ylabel("y")

plt.show()
