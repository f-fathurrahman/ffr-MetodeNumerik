import numpy as np

def f(X):
    x1 = X[0]
    x2 = X[1]
    x3 = X[2]
    f1 = 3*x1 - np.cos(x2*x3) - 1.5
    f2 = 4*x1**2 - 625*x2**2 + 2*x3 - 1
    f3 = 20*x3 + np.exp(-x1*x2) + 9
    return np.array([f1,f2,f3])

def calc_jac(X):
    x1 = X[0]
    x2 = X[1]
    x3 = X[2]
    #
    J11 = 3.0
    J12 = x3*np.sin(x2*x3)
    J13 = x2*np.sin(x2*x3)
    #
    J21 = 8*x1
    J22 = -1250*x2
    J23 = 2.0
    #
    J31 = -x2*np.exp(-x1*x2)
    J32 = -x1*np.exp(-x1*x2)
    J33 = 20.0
    return np.array([
        [J11, J12, J13],
        [J21, J22, J23],
        [J31, J32, J33]
    ])

X = np.array([1.0, 1.0, 1.0])
for i in range(1,50):
    fX = f(X)
    nfX = np.linalg.norm(fX)
    print("X = ", X, "nfX = ", nfX)
    if nfX <= 1e-10:
        print("Converged")
        break
    # Jacobian matrix elements
    J = calc_jac(X)
    invJ = np.linalg.inv(J)
    # Update X
    Xnew = X - np.matmul(invJ, fX)
    X = np.copy(Xnew)
