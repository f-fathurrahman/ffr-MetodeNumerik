import numpy as np

A = np.array([
    [888445.0, 887112.0],
    [887112.0, 885781.0]
])

b = np.array([10.0, 0.0])

x = np.linalg.solve(A, b)
print("x = ", x)

print("residual = ", A @ x - b)

# Manual

x0 = (-A[0,1]*b[1] + A[1,1]*b[0])/(A[0,0]*A[1,1] - A[0,1]*A[1,0])
x1 = (A[0,0]*b[1] - A[1,0]*b[0])/(A[0,0]*A[1,1] - A[0,1]*A[1,0])

print("x0 = ", x0)
print("x1 = ", x1)