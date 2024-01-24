import numpy as np
from LU_decomp import *

A = np.matrix([[6, -4, 1], [-4, 6, -4], [1, -4, 6]], dtype=np.float64)
print(A)

# need to use transpose to create a column matrix
b = np.matrix([-14, 36, 6], dtype=np.float64).transpose()
print(b)

x = LU_solve(LU_decomp(A), b)

print("x = ")
print(x)

print("Check")
print(A*x - b)
