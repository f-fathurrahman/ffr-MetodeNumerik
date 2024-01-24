"""
4x_{1} - x_{2} = 9 \\
-x_{i-1} + 4x_{i} - x_{i+1} = 5,\,\,\, i = 2,\ldots,N-1
-x_{N-1} + 4x_{N} = 5
"""

import numpy as np

N = 50
A = np.matrix(np.zeros((N,N)))

A[0,0] = 4.0
A[0,1] = -1
for i in range(1,N-1):
    A[i,i] = 4.0
    A[i,i-1] = -1.0
    A[i,i+1] = -1.0
A[N-1,N-2] = -1.0
A[N-1,N-1] = 4.0

print(np.linalg.inv(A))
