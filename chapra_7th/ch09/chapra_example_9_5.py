import numpy as np
from linsolve_gauss_elim import *

# Solve Ax = b

A = np.array([
    [3.0, -0.1, -0.2],
    [0.1, 7.0, -0.3],
    [0.3, -0.2, 10.0]
])

b = np.array([
    [7.85],
    [-19.3],
    [71.4]
])

x = linsolve_gauss_elim(A, b)
print("\nSolusi x = ")
print(x)