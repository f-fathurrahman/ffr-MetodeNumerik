import numpy as np
import matplotlib.pyplot as plt

def fit_multilinear_2(x1, x2, y):
    N = len(x1)
    assert N == len(x2)
    assert N == len(y)
    # Build the matrix for linear system
    A = np.zeros((3,3))
    #
    A[0,0] = N
    A[0,1] = np.sum(x1)
    A[0,2] = np.sum(x2)
    #
    A[1,0] = A[0,1]
    A[1,1] = np.sum(x1**2)
    A[1,2] = np.sum(x1*x2)
    #
    A[2,0] = A[0,2]
    A[2,1] = A[1,2]
    A[2,2] = np.sum(x2**2)
    #
    # The RHS vector
    #
    b = np.zeros(3)
    b[0] = np.sum(y)
    b[1] = np.sum(x1*y)
    b[2] = np.sum(x2*y)
    #
    # Solve the linear system
    #
    xsol = np.linalg.solve(A,b)
    return xsol[0], xsol[1], xsol[2]

data = np.loadtxt("table_17_5.dat")
x1 = data[:,0]
x2 = data[:,1]
y  = data[:,2]
a0, a1, a2 = fit_multilinear_2(x1, x2, y)

print("a0 = ", a0)
print("a1 = ", a1)
print("a2 = ", a2)