import numpy as np
from gauss_elim_pivot import *

A = np.matrix([[2.34, -4.10, 1.78], [1.98, 3.47, -2.22], [2.36, -15.17, 6.81]])
print("Matrix A")
print(A)

# need to use transpose to create a column matrix
b = np.matrix([[0.02, -0.73, -6.63]]).transpose()
print("Matrix B")
print(b)

x = gauss_elim_pivot(A, b)

print("x = ")
print(x)

print("Error:")
print(A*x - b)