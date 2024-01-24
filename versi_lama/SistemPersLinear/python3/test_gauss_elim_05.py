import numpy as np
from gauss_elim_pivot import *

A = np.matrix([[4, -2, 1], [-2, 1, -1], [-2, 3, 6]], dtype=np.float64)
print("Matrix A")
print(A)

# need to use transpose to create a column matrix
b = np.matrix([[2, -1, 0]], dtype=np.float64).transpose()
print("Matrix B")
print(b)

x = gauss_elim_pivot(A, b)

print("x = ")
print(x)

print("Error:")
print(A*x - b)