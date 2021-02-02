# Solves Navier-Stokes equation for orifice flow
import numpy as np
import matplotlib.pyplot as plt

Niter = 700
Ndown = 20
Nx = 17
N2x = 2*Nx
Ny = 156
Nb = 15
h = 0.4
h2 = h*h
g = 980.0
nu = 0.5
iiter = 0;
Vtop = 8.0e-4;
omega = 0.1
R = Vtop*h/nu

u = np.zeros( (Nx+1, Ny+1) )
w = np.zeros( (Nx+1, Ny+1) )
ua = np.zeros( (N2x, Ny) )
Torri = open("Torri.dat", "w")
uall = open("uall.dat", "w")


def BelowHole():
    for i in range(Nb+1,Nx+1): # Below orifice
        u[i,0] = u[i-1,1] # du/dy = vx = 0
        w[i-1,0] = w[i-1,1] # Water is at floor
        for j in range(Ndown+1):
            if i == Nb:
                vy = 0
            if i == Nx:
                vy = -np.sqrt(2.0*g*h*(Ny+Nb-j))
            if i == Nx-1:
                vy = -np.sqrt(2.0*g*h*(Ny+Nb-j))/2.0
            u[i,j] = u[i-1,j] - vy*h  # du/ dx= -vy


def BorderRight():
    # Center orifice very sensitive
    for j in range(1,Ny+1):
        vy = -np.sqrt(2.0*g*h*(Ny-j))
        u [Nx,j] = u[Nx-1,j] + vy*h
        u [Nx,j] = u[Nx,j-1]
        w [Nx,j] = -2*( u[Nx,j] - u[Nx,j-1] ) / h**2


def BottomBefore():
    for i in range(1,Nb+1): # Bottom, before the hole
        u[i,Ndown] = u[i,Ndown-1]
        w[i,Ndown] = -2*( u[i,0] - u[i,1] )

def Top():
    for i in range(1,Nx):
        u[i,Ny] = u[i,Ny-1]
        w[i,Ny] = 0.0


def Left () :
    for j in range(Ndown,Ny): # Left wall
        w[0,j] = -2*( u[0,j] - u[1,j] ) / h**2
        u[0,j] = u[1,j]

# du/dx = 0
def Borders(iiter):
    BelowHole()
    BorderRight()
    BottomBefore()
    Top()
    Left()

def Relax(iiter):
    Borders(iiter)
    for i in range(1,Nx):
        for j in range(1,Ny):
            if j <= Ndown:
                if i > Nb:
                    r1 = omega*( ( u[i+1,j] + u[i-1,j] + u[i,j+1] + u[i,j-1] + h*h*w[i,j] )/4 - u[i,j] )
                    u[i,j] += r1
            if j > Ndown:
                r1 = omega * ( ( u [i+1,j] + u[i-1,j] + u[i,j+1] + u[i,j-1] + h*h*w[i,j] )/4 - u[i,j] )
                u[i,j] += r1
    if iiter%50 == 0:
        print("Residual r1: ", r1)
    
    Borders(iiter)
    
    # Relax stream function
    for i in range(1,Nx): 
        for j in range(1,Ny):
            if j <= Ndown:
                if i >= Nb:
                    a1 = w[i+1,j] + w[i-1,j] + w[i,j+1] + w[i,j-1]
                    a2 = ( u[i,j+1] - u[i,j-1] ) * (w[i +1 , j ] - w[i-1,j] )
                    a3 = ( u[i+1,j] - u[i-1,j] ) * (w[i , j +1] - w[i,j-1] )
                    r2 = omega * ( ( a1 + (R/4.0) * ( a3 - a2 ) ) / 4.0 - w[i,j] )
                    w[i,j] += r2
            if j > Ndown:
                a1 = w [i+1,j] + w[i-1,j] + w[i,j+1] + w[i,j-1]
                a2 = ( u[i,j+1] - u[i,j-1] ) * ( w[i+1,j] - w[i-1,j] )
                a3 = ( u[i+1,j] - u[i-1,j] ) * ( w[i,j+1] - w[i,j-1] )
                r2 = omega * ( ( a1 + (R/4.0) * ( a3 - a2 ) ) /4.0 - w[i,j] )
                w[i,j] += r2

while(iiter <= Niter):
    if iiter%100 == 0:
        print("Iteration: ", iiter)
    Relax(iiter)
    iiter += 1


# Send w to disk in gnuplot format
#for j in range(Ny):
#    for i in range(Nx):
#        Torri.write("%8.3e\n" % (w[i,j]))
#    Torri.write("\n")
#Torri.close()

# Send symmetric tank data to disk
#for j in range(Ny): 
#    for i in range(N2x):
#        if i <= Nx:
#            ua[i,j] = u[i,j]
#            uall.write("%8.3e\n" % (ua[i,j]))
#        if i > Nx:
#            ua[i,j] = u[N2x-i,j]
#            uall.write("%8.3e\n" % (ua[i,j]))
#    uall.write("\n")
#uall.close()

#utorr = open("Torri.dat", "w")
## Send u data to disk
#for j in range(Ny):
#    utorr.write("\n")
#    for i in range(Nx):
#        utorr.write("%10.3e\n" % (u[i,j]))
#utorr.close()

x = np.arange(Nx-1)
y = np.arange(Ny-1)
X, Y = np.meshgrid(x,y)
Z = u[X,Y]
plt.clf()
plt.contour(X, Y, Z, levels=20)
plt.savefig("IMG_torri_u_v1.pdf")

x = np.arange(Nx-1)
y = np.arange(Ny-1)
X, Y = np.meshgrid(x,y)
Z = w[X,Y]
plt.clf()
plt.contour(X, Y, Z, levels=20)
plt.savefig("IMG_torri_w_v1.pdf")

