import numpy as np
import matplotlib.pyplot as plt

def u_exact(x, y):
    return np.exp(-(x**2 + y**2))*np.sin(2*x*y)

def f_func(x, y):
    return (-16*x*y*np.cos(2*x*y) + np.sin(2*x*y)*np.cos(2*y) - 4*np.sin(2*x*y))*np.exp(-x**2 - y**2)

def g_func(x, y):
    return np.cos(2*y)

def bx0(y, x0):
    return u_exact(x0, y)

def bxf(y, xf):
    return u_exact(xf, y)

def by0(x, y0):
    return u_exact(x, y0)

def byf(x, yf):
    return u_exact(x, yf)

def calc_sum_of_bv(u, Nx, Ny):
    return np.sum(u[0,:]) + np.sum(u[Nx,:]) + np.sum(u[:,0]) + np.sum(u[:,Ny])

x0 = -np.pi/2; xf = np.pi/2
y0 = -np.pi/2; yf = np.pi/2

Nx = 75; Ny = 75
TOL = 1e-6
NMaxIter = 10000

dx = (xf - x0)/Nx
x = np.linspace(x0, xf, Nx+1)

dy = (yf - y0)/Ny
y = np.linspace(y0, yf, Ny+1)

X, Y = np.meshgrid(x, y)

Nx1 = Nx + 1
Ny1 = Ny + 1

u = np.zeros( (Nx+1,Ny+1) )

for j in range(Ny1):
    u[0,j] = bx0(y[j], x0)
    u[Nx,j] = bxf(y[j], xf)

for i in range(Nx1):
    u[i,0] = by0(x[i], y0)
    u[i,Ny] = byf(x[i], yf)


sum_of_bv = calc_sum_of_bv(u, Nx, Ny)
print("sum_of_bv = ", sum_of_bv)
u[1:Nx,1:Ny] = sum_of_bv/(2*(Nx + Ny - 2))

dx2 = dx**2; dy2 = dy**2; dxy2 = 2*(dx2 + dy2)
rx = dx2/dxy2; ry = dy2/dxy2; rxy = dx2*dy2/dxy2

g = g_func(X, Y)
f = f_func(X, Y)

u0 = np.copy(u)
is_converged = False
for iiter in range(NMaxIter):
    for i in range(1,Nx):
        for j in range(1,Ny):
            u[i,j] = rx*( u[i,j+1] + u[i,j-1] ) + \
                     ry*( u[i+1,j] + u[i-1,j] ) + \
                     rxy*( g[i,j]*u[i,j] - f[i,j] )
    max_error = np.max(np.abs(u0 - u))
    print("Iterasi = %d, max_error = %18.10e" % (iiter, max_error))
    if max_error < TOL:
        is_converged = True
        print("Konvergen pada iterasi = %d, max_error = %18.10e" % (iiter, max_error))
        break
    u0 = np.copy(u)

sum_of_bv = calc_sum_of_bv(u, Nx, Ny)
print("sum_of_bv after iteration = ", sum_of_bv)

if not is_converged:
    print("WARNING: Iterasi tidak konvergen dalam NMaxIter = ", NMaxIter)

fig = plt.figure()
ax = fig.add_subplot(121, projection="3d")
ax2 = fig.add_subplot(122, projection="3d")

du = u_exact(X, Y) - u.T
du_norm = np.dot(du.flatten(), du.flatten())/len(du)
print("du_norm = ", du_norm)

ax.plot_surface(X, Y, u.T, cmap="coolwarm")
ax.set_title("Solusi numerik")
ax.set_xlabel("x")
ax.set_ylabel("y")

ax2.plot_surface(X, Y, u_exact(X,Y), cmap="coolwarm")
ax2.set_title("Solusi analitik")
ax2.set_xlabel("x")
ax2.set_ylabel("y")


plt.show()
