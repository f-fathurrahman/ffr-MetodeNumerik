import numpy as np
import matplotlib.pyplot as plt

def u_exact(x, y):
    return 1/np.sinh(2*np.pi) * np.sin(2*x) * np.sinh(2*(np.pi-y))

bx0 = lambda y: 0.0
bxf = lambda y: 0.0
by0 = lambda x: np.sin(2*x)
byf = lambda x: 0.0

x0 = 0.0; xf = np.pi
y0 = 0.0; yf = np.pi

Nx = 75; Ny = 75
TOL = 1e-6
NMaxIter = 10000

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


sum_of_bv = np.sum(u[0,:]) + np.sum(u[Nx,:]) + \
            np.sum(u[:,0]) + np.sum(u[:,Ny])
u[1:Nx,1:Ny] = sum_of_bv/(2*(Nx + Ny - 2))
#print(np.sum(abs(u[1:Nx,1:Ny]))) # should be close to zero
#print(np.sum(abs(u))) # should not be zero


dx2 = dx*dx; dy2 = dy*dy; dxy2 = 2*(dx2 + dy2)
rx = dx2/dxy2; ry = dy2/dxy2; rxy = rx*dy2

g = np.zeros( (Nx+1,Ny+1) )
f = np.zeros( (Nx+1,Ny+1) )

u0 = np.copy(u)
is_converged = False
for iiter in range(NMaxIter):
    for i in range(1,Nx):
        for j in range(1,Ny):
            u[i,j] = ry*( u[i,j+1] + u[i,j-1] ) + \
                     rx*( u[i+1,j] + u[i-1,j] ) + \
                     rxy*( g[i,j]*u[i,j] - f[i,j] )
    max_error = np.max(abs(u0 - u))
    #print("Iterasi = %d, max_error = %18.10e" % (iiter, max_error))
    if max_error < TOL:
        is_converged = True
        print("Konvergen pada iterasi = %d, max_error = %18.10e" % (iiter, max_error))
        break
    u0 = np.copy(u)

if not is_converged:
    print("WARNING: Iterasi tidak konvergen dalam NMaxIter = ", NMaxIter)

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
X, Y = np.meshgrid(x, y)

du = u_exact(X, Y) - u.T
du_norm = np.dot(du.flatten(), du.flatten())/len(du)
print("du_norm = ", du_norm)

# Notice the transpose for u
ax.plot_surface(X, Y, u.T, cmap="coolwarm")

ax.set_title("Solusi numerik")
ax.set_xlabel("x")
ax.set_ylabel("y");
plt.show()
