import numpy as np
from gauss_elim_pivot import *

A = np.matrix([
    [1.3174, 2.7250, 2.7250, 1.7181],
    [0.4002, 0.8278, 1.2272, 2.5322],
    [0.8218, 1.5608, 0.3629, 2.9210],
    [1.9664, 2.0011, 0.6532, 1.9945]])

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
