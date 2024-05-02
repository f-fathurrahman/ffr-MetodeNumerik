# Solves the 2D Poisson equation in a rectangular domain
from math import *

from numba import jit
import numba
numba.config.DISABLE_JIT = True

@jit(nopython=True)
def PoissonXY(u, x, y, nx, ny, Func, CondX, CondY, itmax=10000, TOL=1e-5):

    f = [[0]*(ny+1) for i in range(nx+1)]
    betXmin = [0]*(ny+1); betXmax = [0]*(ny+1)
    gamXmin = [0]*(ny+1); gamXmax = [0]*(ny+1)
    betYmin = [0]*(nx+1); betYmax = [0]*(nx+1)
    gamYmin = [0]*(nx+1); gamYmax = [0]*(nx+1)

    hx = (x[nx]-x[1])/(nx-1); kx = 1e0/(hx*hx)                 # mesh spacings
    hy = (y[ny]-y[1])/(ny-1); ky = 1e0/(hy*hy) 
    kxy = 2e0*(kx + ky)

    for j in range(2,ny):
        for i in range(2,nx):
            f[i][j] = Func(x[i],y[j])

    for i in range(1,nx+1):
        (alf_min,bet_min,gam_min,alf_max,bet_max,gam_max) = CondY(x[i])
        betYmin[i] = bet_min/(alf_min*hy + bet_min)
        gamYmin[i] = gam_min/(alf_min + bet_min/hy)
        betYmax[i] = bet_max/(alf_max*hy + bet_max)
        gamYmax[i] = gam_max/(alf_max + bet_max/hy)

    for j in range(2,ny):
        (alf_min,bet_min,gam_min,alf_max,bet_max,gam_max) = CondX(y[j])
        betXmin[j] = bet_min/(alf_min*hx + bet_min)
        gamXmin[j] = gam_min/(alf_min + bet_min/hx)
        betXmax[j] = bet_max/(alf_max*hx + bet_max)
        gamXmax[j] = gam_max/(alf_max + bet_max/hx)

    err = 0.0
    for it in range(1,itmax+1):                  # Gauss-Seidel iteration loop
        j = 1                                                  # lower boundary
        for i in range(1,nx+1):
            uij = betYmin[i]*u[i][2] + gamYmin[i]
            eij = 1e0 - u[i][j]/uij if uij else uij - u[i][j]
            if (fabs(eij) > err): err = fabs(eij)
            u[i][j] = uij

        for j in range(2,ny):
            i = 1                                                # left boundary
            uij = betXmin[j]*u[i+1][j] + gamXmin[j]
            eij = 1e0 - u[i][j]/uij if uij else uij - u[i][j]
            if (fabs(eij) > err):
                err = fabs(eij)
            u[i][j] = uij

            for i in range(2,nx):                         # interior mesh points
                uij = (kx*(u[i-1][j] + u[i+1][j]) + ky*(u[i][j-1] + u[i][j+1]) - f[i][j]) / kxy
                eij = 1e0 - u[i][j]/uij if uij else uij - u[i][j]   # local error
                if fabs(eij) > err:
                    err = fabs(eij)             # maximum error
                u[i][j] = uij

            i = nx                                              # right boundary
            uij = betXmax[j]*u[i-1][j] + gamXmax[j]
            eij = 1e0 - u[i][j]/uij if uij else uij - u[i][j]
            if (fabs(eij) > err):
                err = fabs(eij)
            u[i][j] = uij

        j = ny                                                 # upper boundary
        for i in range(1,nx+1):
            uij = betYmax[i]*u[i][ny-1] + gamYmax[i]
            eij = 1e0 - u[i][j]/uij if uij else uij - u[i][j]
            if (fabs(eij) > err):
                err = fabs(eij)
            u[i][j] = uij

        if err <= eps:
            print("Converged")
            break

    if it >= itmax:
        print("PoissonXY: max. number of iterations exceeded !")

    print("err = ", err)
    print("it = ", it)
    return 0


@jit(nopython=True)
def Func(x, y):
    return -(cos(x+y) + cos(x-y))

@jit(nopython=True)
def CondX(y):
    alf_min = 1e0; bet_min = 0e0; gam_min = -cos(y)
    alf_max = 1e0; bet_max = 0e0; gam_max = -cos(y)
    return (alf_min, bet_min, gam_min, alf_max, bet_max, gam_max)

@jit(nopython=True)
def CondY(x):
    alf_min = 1e0; bet_min = 0e0; gam_min = 0e0
    alf_max = 1e0; bet_max = 0e0; gam_max = 0e0
    return (alf_min, bet_min, gam_min, alf_max, bet_max, gam_max)

# main

xmin = -pi; xmax = pi; ymin = -pi/2; ymax = pi/2          # domain boundaries
nx = 51; ny = 26                                      # number of mesh points
eps = 1e-5                                      # relative solution tolerance

u = [[0]*(ny+1) for i in range(nx+1)]                              # solution
x = [0]*(nx+1); y = [0]*(ny+1)                       # mesh point coordinates

hx = (xmax-xmin)/(nx-1)
for i in range(1,nx+1): x[i] = xmin + (i-1)*hx                # x-mesh points
hy = (ymax-ymin)/(ny-1)
for j in range(1,ny+1): y[j] = ymin + (j-1)*hy                # y-mesh points

for j in range(1,ny+1):               # initial approximation of the solution
    for i in range(1,nx+1): u[i][j] = 0e0

import time
start = time.perf_counter()
PoissonXY(u, x, y, nx, ny, Func, CondX, CondY, TOL=eps)
end = time.perf_counter()
print(f"Elapsed time = {end - start} s")

