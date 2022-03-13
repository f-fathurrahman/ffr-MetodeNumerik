import numpy as np
from linsolve_gauss_elim_verbose import *

# Solve Ax = b

A = np.array([
    [70, 1.0, 0.0],
    [60, -1.0, 1.0],
    [40.0, 0.0, -1.0]
])

b = np.array([
    [636.7],
    [518.6],
    [307.4]
])

x = linsolve_gauss_elim(A, b, verbose=True)
print("Solution x = ")
print(x)

print("Check: (should be close to zero)")
print(A @ x - b)