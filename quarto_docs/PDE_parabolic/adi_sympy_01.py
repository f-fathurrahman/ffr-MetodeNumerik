from sympy import *

Nx = 4
Ny = 5

uk = MatrixSymbol("u^{k}", Nx+1, Ny+1)
ukp1 = MatrixSymbol("u^{k+1}", Nx+1, Ny+1)

