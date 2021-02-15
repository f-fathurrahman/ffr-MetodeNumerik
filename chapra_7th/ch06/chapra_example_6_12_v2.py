import numpy as np

def f(X):
    x = X[0]
    y = X[1]
    f1 = x**2 + x*y - 10
    f2 = y + 3*x*y**2 - 57
    return np.array([f1,f2])

def calc_jac(X):
    x = X[0]
    y = X[1]
    dudx = 2*x + y
    dudy = x
    dvdx = 3*y**2
    dvdy = 1 + 6*x*y
    return np.array([
        [dudx, dudy],
        [dvdx, dvdy]
    ])

X = np.array([1.5, 3.5])
for i in range(1,6):
    fX = f(X)
    nfX = np.linalg.norm(fX)
    print("X = ", X, "nfX = ", fX)
    if nfX <= 1e-10:
        print("Converged")
        break
    # Jacobian matrix elements
    J = calc_jac(X)
    invJ = np.linalg.inv(J)
    # Update X
    Xnew = X - np.matmul(invJ, fX)
    X = np.copy(Xnew)
