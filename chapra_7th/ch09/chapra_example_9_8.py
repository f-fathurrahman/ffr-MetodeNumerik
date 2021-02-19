import numpy as np
from linsolve_gauss_elim_v2 import *

# Solve Ax = b

A = np.array([
    [1.0, 0.667],
    [-0.5, 1],
])

b = np.array([
    [6.0],
    [1.0],
])

linsolve_gauss_elim(A, b)
#x = linsolve_gauss_elim(A, b)
#print("\nSolusi x = ")

#print(x)