import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

Nmax = 100
Niter = 600
V = np.zeros((Nmax,Nmax))
grid = np.ones((Nmax,Nmax))
sq3 = np.sqrt(3.0)
y0 = -sq3/4.0

def contour():
    for j in range(Nmax):
        V[0,j] = 1000.0


# Set interior grid points to 0
for i in range(Nmax):
    y = y0 + i*0.01
    x0 = -0.5
    for j in range(Nmax):
        x = x0 + j*0.01
        cond1 = (y <= sq3*(x+0.25) and x < 0.0) or (y < -sq3*(x - 0.25) and x >= 0)
        if cond1:
            grid[i,j] = 0.0
        else:
            if(y <= sq3/4.0 + 0.01):
                V[i,j] = 0.0 # Triangletip


for iiter in range(1,Niter):
    if(iiter%50==0):
        print("Doing iteration = ", iiter)
    contour() # keep one size at 1000 V (apply boundary condition)
    for i in range(1,Nmax-2):
        for j in range(1,Nmax-2):
            # XXX FIXME XXX
            if grid[i,j] == 0.0:
                V[i,j] = 0.25*(V[i+1,j] + V[i-1,j] + V[i,j+1] + V[i,j-1])

x = range(0,Nmax-1,2)
y = range(0,Nmax,2)
X, Y = np.meshgrid(x,y)

Z = V[X,Y]

fig = plt.figure()
ax = Axes3D(fig)
ax.plot_wireframe(X, Y, Z, color='b')
ax.set_xlabel("X")
ax.set_xlabel("Y")
ax.set_xlabel("Potential")
plt.savefig("IMG_LaplaceTri.pdf")
