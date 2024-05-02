# Solves the 2D Poisson equation in a rectangular domain
import math
import numpy as np

from numba import jit
import numba
numba.config.DISABLE_JIT = False

@jit(nopython=True)
def PoissonXY(u, x, y, Nx, Ny, rhs_func, bcond_x, bcond_y, NiterMax=10000, TOL=1e-5):

    f = np.zeros((Nx+1,Ny+1))

    betXmin = np.zeros(Ny+1)
    betXmax = np.zeros(Ny+1)
    gamXmin = np.zeros(Ny+1)
    gamXmax = np.zeros(Ny+1)
    
    betYmin = np.zeros(Nx+1)
    betYmax = np.zeros(Nx+1)
    gamYmin = np.zeros(Nx+1)
    gamYmax = np.zeros(Nx+1)

    hx = (x[Nx]-x[1])/(Nx-1)
    kx = 1.0/(hx*hx)
    
    hy = (y[Ny]-y[1])/(Ny-1)
    ky = 1.0/(hy*hy) 
    kxy = 2.0*(kx + ky)

    # For interior points
    for j in range(2,Ny):
        for i in range(2,Nx):
            f[i,j] = rhs_func(x[i], y[j])

    for i in range(1,Nx+1):
        alf_min, bet_min, gam_min, alf_max, bet_max, gam_max = bcond_y(x[i])
        betYmin[i] = bet_min/(alf_min*hy + bet_min)
        gamYmin[i] = gam_min/(alf_min + bet_min/hy)
        betYmax[i] = bet_max/(alf_max*hy + bet_max)
        gamYmax[i] = gam_max/(alf_max + bet_max/hy)

    for j in range(2,Ny):
        alf_min, bet_min, gam_min, alf_max, bet_max, gam_max = bcond_x(y[j])
        betXmin[j] = bet_min/(alf_min*hx + bet_min)
        gamXmin[j] = gam_min/(alf_min + bet_min/hx)
        betXmax[j] = bet_max/(alf_max*hx + bet_max)
        gamXmax[j] = gam_max/(alf_max + bet_max/hx)

    iterNo = 1
    err = 0.0

    while iterNo <= NiterMax:
        # lower boundary
        j = 1
        for i in range(1,Nx+1):
            uij = betYmin[i]*u[i,2] + gamYmin[i]
            eij = 1e0 - u[i,j]/uij if uij else uij - u[i,j]
            if math.fabs(eij) > err:
                err = math.fabs(eij)
            u[i,j] = uij

        for j in range(2,Ny):
            # left boundary
            i = 1
            uij = betXmin[j]*u[i+1,j] + gamXmin[j]
            eij = 1.0 - u[i,j]/uij if uij else uij - u[i,j]
            if math.fabs(eij) > err:
                err = math.fabs(eij)
            u[i,j] = uij

            # interior mesh points
            for i in range(2,Nx):
                uij = ( kx*(u[i-1,j] + u[i+1,j]) + ky*(u[i,j-1] + u[i,j+1]) - f[i,j]) / kxy
                eij = 1.0 - u[i,j]/uij if uij else uij - u[i,j]
                if math.fabs(eij) > err:
                    err = math.fabs(eij)
                u[i,j] = uij

            # right boundary
            i = Nx
            uij = betXmax[j]*u[i-1,j] + gamXmax[j]
            eij = 1.0 - u[i,j]/uij if uij else uij - u[i,j]
            if math.fabs(eij) > err:
                err = math.fabs(eij)
            u[i,j] = uij

        # upper boundary
        j = Ny
        for i in range(1,Nx+1):
            uij = betYmax[i]*u[i,Ny-1] + gamYmax[i]
            eij = 1.0 - u[i,j]/uij if uij else uij - u[i,j]
            if math.fabs(eij) > err:
                err = math.fabs(eij)
            u[i,j] = uij
        #
        if err <= TOL:
            print("Converged")
            break
        #
        iterNo += 1

    if iterNo >= NiterMax:
        print("WARNING: maximum number of iterations exceeded !")
    
    return iterNo, err


@jit(nopython=True)
def rhs_func(x, y):
    return -math.cos(x+y) - math.cos(x-y)


# Coefficients for left and right boundaries
@jit(nopython=True)
def bcond_x(y):
    alf_min = 1.0; bet_min = 0.0; gam_min = -math.cos(y)
    alf_max = 1.0; bet_max = 0.0; gam_max = -math.cos(y)
    return alf_min, bet_min, gam_min, alf_max, bet_max, gam_max

# Coefficients for lower and upper boundaries
@jit(nopython=True)
def bcond_y(x):
    alf_min = 1.0; bet_min = 0.0; gam_min = 0.0
    alf_max = 1.0; bet_max = 0.0; gam_max = 0.0
    return alf_min, bet_min, gam_min, alf_max, bet_max, gam_max

# main

xmin = -math.pi
xmax = math.pi
ymin = -math.pi/2
ymax = math.pi/2
Nx = 51
Ny = 26
TOL = 1e-5

u = np.zeros((Nx+1, Ny+1))
x = np.zeros(Nx+1)
y = np.zeros(Ny+1)

hx = (xmax-xmin)/(Nx-1)
for i in range(1,Nx+1):
    x[i] = xmin + (i-1)*hx

hy = (ymax-ymin)/(Ny-1)
for j in range(1,Ny+1):
    y[j] = ymin + (j-1)*hy


import time
start = time.perf_counter()
iterNo, err = PoissonXY(u, x, y, Nx, Ny, rhs_func, bcond_x, bcond_y, TOL=TOL, NiterMax=10_000)
end = time.perf_counter()
print(f"iterNo = {iterNo} err={err}")
print(f"Elapsed time = {end-start} s")


import matplotlib.pyplot as plt
XX, YY = np.meshgrid(x[1:Nx+1], y[1:Ny+1])
uu = u[1:Nx+1,1:Ny+1]
plt.clf()
plt.contourf(XX, YY, uu.T)
plt.gca().set_aspect("equal")
plt.show()
