# Solves the 2D Poisson equation in a rectangular domain
import math
import numpy as np

from numba import jit
import numba
numba.config.DISABLE_JIT = False


@jit(nopython=True)
def Poisson0(u, x, y, Nx, Ny, rhs_func, NiterMax=10000, TOL=1e-5):
#----------------------------------------------------------------------------
#  Solves the 2D Poisson equation in Cartesian coordinates with Dirichlet
#  boundary conditions on a regular grid with (nx x ny) nodes (x[],y[]) using
#  the Gauss-Seidel method. The solution u[][] is converged with relative
#  precision eps. An error index is returned: 0 - normal execution.
#  Calls: Func(x,y) - RHS of Poisson equation
#----------------------------------------------------------------------------

    SMALL = 2.220446049250313e-16
    # SMALL = np.finfo(np.float64).eps

    hx = (x[Nx]-x[1])/(Nx-1)
    kx = 1/(hx*hx)
    
    hy = (y[Ny]-y[1])/(Ny-1)
    ky = 1/(hy*hy) 
    
    kxy = 2*(kx + ky)

    # RHS of PDE for interior points
    f = np.zeros((Nx+1,Ny+1))
    for j in range(2,Ny):
        for i in range(2,Nx):
            f[i,j] = rhs_func(x[i], y[j])

    iterNo = 1
    while iterNo <= NiterMax:
        err = 0.0
        for j in range(2,Ny):
            for i in range(2,Nx): # interior mesh points
                uij = (kx*(u[i-1,j] + u[i+1,j]) + ky*(u[i,j-1] + u[i,j+1]) - f[i,j]) / kxy
                if abs(uij) > SMALL:
                    eij = 1.0 - u[i,j]/uij
                else:
                    eij = uij - u[i,j]   # local error
                # Take maximum error
                if math.fabs(eij) > err:
                    err = math.fabs(eij)
                u[i,j] = uij

        if err <= TOL:
            print("Converged")
            # @njit or @jit(nopython=True) does not support formatted print (?)
            break
        
        iterNo += 1

    if iterNo > NiterMax:
        print("WARNING: Number of iterations exceed NiterMax")

    return iterNo, err # return some info


# RHS function for Poisson0
@jit(nopython=True)
def rhs_function(x, y):
    return math.cos(x+y) - math.cos(x-y)

# main

# Define domain boundaries
xmin = -math.pi; xmax = math.pi
ymin = -math.pi; ymax = math.pi
# number of mesh points
Nx = 51
Ny = 51
# relative solution tolerance
TOL = 1e-5

# solution
u = np.zeros((Nx+1,Ny+1))
# Mesh points coordinates
x = np.zeros(Nx+1)
y = np.zeros(Ny+1)
# NOTE: We don't use the 0-th index

hx = (xmax-xmin)/(Nx-1)
for i in range(1,Nx+1):
    x[i] = xmin + (i-1)*hx

hy = (ymax-ymin)/(Ny-1)
for j in range(1,Ny+1):
    y[j] = ymin + (j-1)*hy

# initial values of u are set to zeros

import time
start = time.perf_counter()
iterNo, err = Poisson0(u, x, y, Nx, Ny, rhs_function, TOL=TOL)
end = time.perf_counter()
print(f"iterNo = {iterNo} err={err}")
print("Elapsed time = {}s".format((end - start)))


"""
out = open("Poisson.txt","w")                              # open output file
out.write("      x         y          u\n")
for j in range(1,ny+1):
    for i in range(1,nx+1):
        out.write(("{0:10.5f}{1:10.5f}{2:14.5e}\n").format(x[i],y[j],u[i][j]))
out.close()

umin = umax = u[1][1]                   # minimum and maximum of the solution
for j in range(1,ny+1):
    for i in range(1,nx+1):
        if (u[i][j] < umin): umin = u[i][j]
        if (u[i][j] > umax): umax = u[i][j]
"""
