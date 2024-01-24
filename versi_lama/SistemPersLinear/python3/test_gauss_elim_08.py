import numpy as np
from gauss_elim_pivot import *

A = np.matrix([
    [0, 2, 5, -1],
    [2, 1, 3, 0],
    [-2, -1, 3, 1],
    [3, 3, -1, 2]], dtype=np.float64)
print("Matrix A")
print(A)

# need to use transpose to create a column matrix
b = np.matrix([[-3, 3, -2, 5]], dtype=np.float64).transpose()
print("Matrix B")
print(b)

x = gauss_elim_pivot(A, b)

print("x = ")
print(x)

print("Error:")
print(A*x - b)