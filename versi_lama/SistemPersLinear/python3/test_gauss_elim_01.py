import numpy as np
from gauss_elim import *

A = np.matrix([[6, -4, 1], [-4, 6, -4], [1, -4, 6]], dtype=np.float64)
print(A)

# need to use transpose to create a column matrix
b = np.matrix([[-14, 36, 6],[22, -18, 7]], dtype=np.float64).transpose()
print(b)

x = gauss_elim(A, b)

print("x = ")
print(x)

print(A*x - b)