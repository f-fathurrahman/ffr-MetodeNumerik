import numpy as np
from cholesky_decomp import *

A = np.matrix([
    [4, -2, 2],
    [-2, 2, -4],
    [2, -4, 11]], dtype=np.float64)
print(A)

b = np.matrix([1, 3/2, 3], dtype=np.float64).transpose()
print(b)

L = cholesky_decomp(A)
print("L after cholesky_decomp:")
print(L)
print("Check L")
print(L*L.transpose())

x = cholesky_solve(L, b)

print("x = ")
print(x)

print("Check")
print(A*x - b)
