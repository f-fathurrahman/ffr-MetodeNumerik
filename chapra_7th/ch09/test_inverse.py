import numpy as np

A = np.array([
    [10, 3 , -9, 6, 4],
    [2 , -1, 6, 7, 1],
    [3 , 2 , -3, 15, 5],
    [8 , -1, 1, 4, 2],
    [11, 1 , -2, 18, 7]
], dtype=np.float64)

Ainv = np.linalg.inv(A)
print("Ainv = ")
print(Ainv)

print("Check: (should be close to zero)")
print(A @ Ainv - np.eye(A.shape[0]))
