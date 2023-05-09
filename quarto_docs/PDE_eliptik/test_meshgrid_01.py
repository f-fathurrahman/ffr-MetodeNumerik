import numpy as np
import matplotlib.pyplot as plt

x0 = 0.5
y0 = -0.1
def u_exact(x, y):
    arg1 = (x-x0)**2 + (y-y0)**2
    return np.exp(-10*arg1)

Nx = 50
x = np.linspace(-1.0, 1.0, Nx)

Ny = 51
y = np.linspace(-1.0, 1.0, Ny)

X, Y = np.meshgrid(x, y)

fig = plt.figure()

ax1 = fig.add_subplot(121, projection="3d")
ax1.plot_surface(X, Y, u_exact(X,Y), cmap="coolwarm")
ax1.set_title("Using np.array meshgrid directly")
ax1.set_xlabel("x")
ax1.set_ylabel("y")

u = np.zeros((Nx,Ny))
for i in range(Nx):
    for j in range(Ny):
        u[i,j] = u_exact(x[i], y[j])

ax2 = fig.add_subplot(122, projection="3d")
ax2.plot_surface(X, Y, u.T, cmap="coolwarm")
ax2.set_title("Manual")
ax2.set_xlabel("x")
ax2.set_ylabel("y")

plt.show()
