import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

Nmax = 100
Niter = 50
V = np.zeros((Nmax,Nmax))

for k in range(Nmax-1):
    V[0,k] = 100.0

for iiter in range(Niter):
    if iiter%10 == 0:
        print("Iteration = ", iiter)
    for i in range(1,Nmax-2):
        for j in range(1,Nmax-2):
            V[i,j] = 0.25*(V[i+1,j] + V[i-1,j] + V[i,j+1] + V[i,j-1])

x = range(0,50,2)
y = range(0,50,2)
X, Y = np.meshgrid(x,y)

Z = V[X,Y]
fig = plt.figure()
ax = Axes3D(fig)
ax.plot_wireframe(X, Y, Z, color='b')
ax.set_xlabel("X")
ax.set_xlabel("Y")
ax.set_xlabel("Potential")
plt.savefig("IMG_LaplaceLine.pdf")
