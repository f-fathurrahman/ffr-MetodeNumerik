from sympy import *

Nx = 3
Ny = 4

D2x = Matrix([
    [-2, 1, 0],
    [1, -2, 1],
    [0, 1, -2]
])

Ix = Matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]
])

D2y = Matrix([
    [-2, 1, 0, 0],
    [1, -2, 1, 0],
    [0, 1, -2, 1],
    [0, 0, 1, -2]
])

Iy = Matrix([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
])

assert Ix.shape[0] == Ix.shape[1]
assert Iy.shape[0] == Iy.shape[1]

D2xy = (KroneckerProduct(D2x,Iy) + KroneckerProduct(Ix,D2y)).doit()

u = MatrixSymbol("u", Nx, Ny)
uvec = Matrix([0]*Nx*Ny)
ip = 0
for i in range(Nx):
    for j in range(Ny):
        uvec[ip] = u[i,j]
        ip += 1

# Boundary points
bcpoints = []
for i in range(Nx):
    for j in range(Ny):
        if (i == 0) or (i == Nx-1):
            bcpoints.append(u[i,j])
        elif (j == 0) or (j == Ny-1):
            bcpoints.append(u[i,j])

