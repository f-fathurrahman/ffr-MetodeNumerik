import numpy as np
from LU_with_pivot import *

A = np.matrix([
    [2, -2, 6],
    [-2, 4, 3],
    [-1, 8, 4]
], dtype=np.float64)

b = np.array([
    [16],
    [0],
    [-1]
])

L, U, iperm = LU_decomp_pivot(A)

print("L = ", L)
print("U = ", U)

x = LU_solve_pivot(L, U, iperm, b)

print("Solusi x = ")
print(x)

print("Cek solusi: Ax - b = ");
print(A@x - b)