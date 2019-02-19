import numpy as np
from gauss_elim_pivot import *

A = np.matrix([[3, -3, 3], [-3, 5, 1], [3, 1, 5]], dtype=np.float64)
print("Matrix A")
print(A)

# need to use transpose to create a column matrix
b = np.matrix([[9, -7, 12]], dtype=np.float64).transpose()
print("Matrix B")
print(b)

x = gauss_elim_pivot(A, b)

print("x = ")
print(x)

print("Error:")
print(A*x - b)