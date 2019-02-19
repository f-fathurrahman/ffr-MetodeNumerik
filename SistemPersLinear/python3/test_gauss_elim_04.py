import numpy as np
from gauss_elim_pivot import *

A = np.matrix([[4, -8, 20], [8, 13, 16], [20, 16, -91]], dtype=np.float64)
print("Matrix A")
print(A)

# need to use transpose to create a column matrix
b = np.matrix([[24, 18, -119]], dtype=np.float64).transpose()
print("Matrix B")
print(b)

x = gauss_elim_pivot(A, b)

print("x = ")
print(x)

print("Error:")
print(A*x - b)