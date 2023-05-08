import numpy as np

α = 1e-2
def initial_cond(X,Y):
    return np.zeros(X.shape)

def bxyt(x,y,t):
    return np.exp(y)*np.cos(x) - np.exp(x)*np.cos(y)

x0 = 0.0
xf = 4.0
y0 = 0.0
yf = 4.0

t0 = 0
tf = 5000
Nx = 40
Ny = 40
NtimeSteps = 50

dt = tf/NtimeSteps
dx = (xf - x0)/Nx
dy = (yf - y0)/Ny

x = np.linspace(x0, xf, Nx+1)
y = np.linspace(y0, yf, Ny+1)
t = np.linspace(t0, tf, NtimeSteps+1)

X, Y = np.meshgrid(x, y)

u = np.zeros((Nx+1,Ny+1))

# Initial condition
u[:,:] = initial_cond(X, Y)


import matplotlib.pyplot as plt

fig = plt.figure()
ax1 = fig.add_subplot(111, projection="3d")
t = 0.0
str_title = "t = {:10.5f}".format(t)
ax1.plot_surface(X, Y, u.T, cmap="coolwarm", edgecolors="k")
ax1.set_title(str_title)
ax1.set_xlabel("x")
ax1.set_ylabel("y")
plt.show()


rx = α**2 * dt / dx**2
ry = α**2 * dt / dy**2

Ay = np.zeros((Nx-1,Nx-1))
for j in range(Nx-1):
    Ay[j,j] = 1 + 2*ry
    if j > 0:
        Ay[j-1,j] = -ry
        Ay[j,j-1] = -ry

Ax = np.zeros((Ny-1,Ny-1))
for i in range(Ny-1):
    Ax[i,i] = 1 + 2*rx
    if i > 0:
        Ax[i-1,i] = -rx
        Ax[i,i-1] = -rx

bx = np.zeros(Nx-1)
by = np.zeros(Ny-1)

u_prev = np.zeros((Nx+1,Ny+1))
for k in range(1,NtimeSteps+1):

    u_prev[:,:] = u[:,:]
    t = k*dt

    # Apply BC 
    for j in range(Ny):
        u[0,j] = bxyt( x[0], y[j], t )
        u[Nx,j] = bxyt( x[Nx], y[j], t )
    
    for i in range(Nx):
        u[i,0] = bxyt( x[i], y[0], t )
        u[i,Ny] = bxyt( x[i], y[Ny], t )
    
    if k % 2 == 1:
        for j in range(1,Ny):
            #print("\nj = ", j)
            for i in range(1,Nx-1):
                #print("Set i for bx = ", i)
                bx[i] = rx*( u_prev[i+1,j] + u_prev[i-1,j] ) + (1 - 2*rx)*u_prev[i,j]
            #
            bx[0] = ry*u[0,j]
            bx[Nx-2] = ry*u[Nx,j]
            u[1:Nx,j] = np.linalg.solve(Ay, bx)
    
    elif k % 2 == 0:
        for i in range(1,Nx):
            #print("\ni = ", i)
            for j in range(1,Ny-1):
                #print("Set i for by = ", i)
                by[j] = ry*( u_prev[i,j+1] + u_prev[i,j-1] ) + (1 - 2*ry)*u_prev[i,j]
            #
            by[0] = rx*u[i,0]
            by[Ny-2] = rx*u[i,Ny]
            u[i,1:Ny] = np.linalg.solve(Ax, by)


fig = plt.figure()
ax1 = fig.add_subplot(111, projection="3d")
str_title = "t = {:10.5f}".format(t)
ax1.plot_surface(X, Y, u.T, cmap="coolwarm", edgecolors="gray")
ax1.set_title(str_title)
ax1.set_xlabel("x")
ax1.set_ylabel("y")
plt.show()
