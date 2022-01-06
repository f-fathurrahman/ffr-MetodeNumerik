import numpy as np
import matplotlib.pyplot as plt

def fit_quadratic_line(x, y):
    N = len(x)
    assert N == len(y)
    # Build the matrix for linear system
    A = np.zeros((3,3))
    #
    A[0,0] = N
    #
    A[0,1] = np.sum(x)
    A[1,0] = A[0,1]
    #
    A[0,2] = np.sum(x**2)
    A[1,1] = A[0,2]
    A[2,0] = A[0,2]
    #
    A[1,2] = np.sum(x**3)
    A[2,1] = A[1,2]
    #
    A[2,2] = np.sum(x**4)
    #
    # The RHS vector
    #
    b = np.zeros(3)
    b[0] = np.sum(y)
    b[1] = np.sum(x*y)
    b[2] = np.sum(x**2 * y)
    #
    # Solve the linear system
    #
    xsol = np.linalg.solve(A,b)
    return xsol[0], xsol[1], xsol[2]

data = np.loadtxt("table_17_4.dat")
x = data[:,0]
y = data[:,1]
a0, a1, a2 = fit_quadratic_line(x, y)

print("a0 = ", a0)
print("a1 = ", a1)
print("a2 = ", a2)