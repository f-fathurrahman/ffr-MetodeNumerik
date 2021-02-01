# Solves Navies-Stokes equation for flow around beam
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import numpy as np

print("Working, look for the figure after 100 iterations")
Nxmax = 70
Nymax = 20
IL = 10
H = 8
T = 8
h = 1.0

u = np.zeros( (Nxmax+1, Nymax+1) ) # Stream
w = np.zeros( (Nxmax+1, Nymax+1) ) # Vorticity
V0 = 1.0
omega = 0.1
nu = 1.0

iter = 0
R = V0*h/nu

def borders():
    for i in range(Nxmax+1):
        for j in range(Nymax+1):
            w[i,j] = 0.0
            u[i,j] = j*V0

    # Fluid surface
    for i in range(Nxmax+1):
        u[i,Nymax] = u[i,Nymax-1] + V0*h
        w[i,Nymax-1] = 0.0

    for j in range(Nymax+1):
        u[1,j] = u[0,j]
        w[0,j] = 0.0
    
    for i in range(Nxmax+1):
        if i <= IL and i >= IL+T:
            u[i,0] = 0.0
            w[i,0] = 0.0
    
    for j in range(1,Nymax):
        w[Nxmax,j] = w[Nxmax-1,j]
        u[Nxmax,j] = u[Nxmax-1,j]


def beam():
    # BC for beam
    for j in range(H+1): # Sides
        w[IL,j] = -2*u[IL-1,j]/(h*h) # Front
        w[IL+T,j] = -2*u[IL+T+1,j]/(h*h) # Back
    #
    for i in range(IL, IL+T+1):
        w[i,H-1] = -2*u[i,H]/(h*h)
    #
    for i in range(IL,IL+T+1):
        for j in range(H+1):
            u[IL,j] = 0.0 # Front
            u[IL+T,j] = 0.0 # Back
            u[i,H] = 0 # top

def relax():
    beam()    # Reset
    for i in range(1,Nxmax): # Relax stream
        for j in range(1,Nymax):
            r1 = omega * ( ( u[i+1,j] + u[i-1,j] + u[i,j+1] + u[i,j-1] + h*h*w[i,j] )/4 - u[i,j] )
            u[i,j] += r1
    #
    for i in range(1,Nxmax): # Relax vorticity
        for j in range(1,Nymax):
            a1 = w [i+1,j] + w[i-1,j] + w[i,j +1] + w[i , j -1]
            a2 = ( u [i,j+1] - u[i,j-1] ) * (w[i+1,j] - w[i-1,j] )
            a3 = ( u [i+1,j] - u[i-1,j] ) * (w[i,j+1] - w[i,j-1] )
            r2 = omega * ( ( a1 - (R/4.0) * ( a2 - a3 ) )/4.0 - w[i,j] )
            w[i,j] += r2


borders() 
iiter = 0
while(iiter <= 100):
    iiter += 1
    if iiter %10 == 0:
        print(iiter)
    relax()

for i in range(0,Nxmax+1):
    for j in range(0,Nymax+1) :
        u[i,j] = u[i,j]/V0/h
x = np.arange(Nxmax-1)
y = np.arange(Nymax-1)
X, Y = np.meshgrid(x,y)

# Stream flow
def functz(u):
    z = u[X,Y]
    return z

Z = u[X,Y]
fig = plt.figure()
ax = Axes3D(fig)
ax.plot_wireframe(X, Y, Z, color="r")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("StreamFunction")
plt.savefig("IMG_beam_v1.pdf")

