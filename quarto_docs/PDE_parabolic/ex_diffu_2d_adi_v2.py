import numpy as np

α = 1.0

def u_exact(x, y, t):
    return np.exp(-t)*np.exp(-x**2 - y**2)

def f_source(x, y, t):
    return (4*x**2 + 4*y**2 - 1)*np.exp(-t)*np.cos(2*x*y)

def initial_cond(x,y):
    return u_exact(x,y,0)

def bx0(y, t, x0):
    return u_exact(x0, y, t)

def bxf(y, t, xf):
    return u_exact(xf, y, t)

def by0(x, t, y0):
    return u_exact(x, y0, t)

def byf(x, t, yf):
    return u_exact(x, yf, t)

x0 = -2.0
xf = 2.0
y0 = -2.0
yf = 2.0

t0 = 0.0
tf = 1.0
Nx = 50
Ny = 50
NtimeSteps = 100

#
dt = tf/NtimeSteps
dx = (xf - x0)/Nx
dy = (yf - y0)/Ny

x = np.linspace(x0, xf, Nx+1)
y = np.linspace(y0, yf, Ny+1)
t = np.linspace(t0, tf, NtimeSteps+1)

X, Y = np.meshgrid(x, y)

u = np.zeros((Nx+1,Ny+1))

# Initial condition
u[:,:] = initial_cond(X, Y).T


import matplotlib.pyplot as plt

fig = plt.figure()
ax1 = fig.add_subplot(111, projection="3d")
t = 0.0
str_title = "t = {:10.5f}".format(t)
ax1.plot_surface(X, Y, u.T, cmap="coolwarm")
ax1.set_title(str_title)
#ax1.set_zlim(0.0, 0.9)
ax1.set_xlabel("x")
ax1.set_ylabel("y")
plt.show()


rx = α * dt / dx**2
rx1 = 1 + 2*rx
rx2 = 1 - 2*rx

ry = α * dt / dy**2
ry1 = 1 + 2*ry
ry2 = 1 - 2*ry

dxyt = dt*dx**2*dy**2

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
    print("k = ", k)

    # Apply BC 
    for j in range(Ny):
        u[0,j] = bx0( y[j], t, x[0] )
        u[Nx,j] = bxf( y[j], t, x[Nx] )
    
    for i in range(Nx):
        u[i,0] = by0( x[i], t, y[0] )
        u[i,Ny] = byf( x[i], t, y[Ny] )
    
    if k % 2 == 1:
        for j in range(1,Ny):
            #print("\nj = ", j)
            for i in range(1,Nx-1):
                #print("Set i for bx = ", i)
                bx[i-1] = rx*( u_prev[i+1,j] + u_prev[i-1,j] ) + (1 - 2*rx)*u_prev[i,j] + dxyt*f_source(x[i], y[j], t)
            #
            bx[0] = ry*u[0,j] - f_source(x[i], y[j], t)*dxyt
            bx[Nx-2] = ry*u[Nx,j] - f_source(x[i], y[j], t)*dxyt
            u[1:Nx,j] = np.linalg.solve(Ay, bx)
    
    elif k % 2 == 0:
        for i in range(1,Nx):
            #print("\ni = ", i)
            for j in range(1,Ny-1):
                #print("Set i for by = ", i)
                by[j-1] = ry*( u_prev[i,j+1] + u_prev[i,j-1] ) + (1 - 2*ry)*u_prev[i,j] + dxyt*f_source(x[i], y[i], t)
            #
            by[0] = rx*u[i,0] - f_source(x[i], y[i], t)*dxyt
            by[Ny-2] = rx*u[i,Ny] - f_source(x[i], y[i], t)*dxyt
            u[i,1:Ny] = np.linalg.solve(Ax, by)

fig = plt.figure()
ax1 = fig.add_subplot(111, projection="3d")
str_title = "t = {:10.5f}".format(t)
ax1.plot_surface(X, Y, u.T, cmap="coolwarm")
ax1.set_title(str_title)
ax1.set_xlabel("x")
ax1.set_ylabel("y")
plt.show()

du = np.abs(u_exact(X, Y, t) - u.T)
print("MAE u = ", np.sum(du)/len(du))
