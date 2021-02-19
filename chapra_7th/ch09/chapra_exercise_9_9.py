import numpy as np
from linsolve_gauss_elim_v2 import *

# Solve Ax = b

A = np.array([
    [8.0, 2.0, -2.0],
    [10.0, 2.0, 4.0],
    [12.0, 2.0, 2.0]
])

b = np.array([
    [-2.0],
    [4.0],
    [6.0]
])

x = linsolve_gauss_elim(A, b)

print(np.matmul(A,x))
#x = linsolve_gauss_elim(A, b)
#print("\nSolusi x = ")

#print(x)