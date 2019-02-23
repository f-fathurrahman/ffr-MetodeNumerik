import numpy as np
from gauss_elim import *

A = np.matrix([[1, 1, 1], [2, 3, 1], [1, -1, -1]], dtype=np.float64)
print(A)

# need to use transpose to create a column matrix
b = np.matrix([4,9,-2], dtype=np.float64).transpose()
print(b)

x = gauss_elim(A, b)

print("x = ")
print(x)

print(A*x - b)