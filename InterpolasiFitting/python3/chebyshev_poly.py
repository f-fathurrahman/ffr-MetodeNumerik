import numpy as np

def transform_to_m1p1(x):
    """
    Pemetaan linear ke [-1,1]
    """
    a = np.min(x)
    b = np.max(x)
    z = (2*x - (b + a))/(b - a)
    return z

def chebyshev_poly_eval(degree, z):
    """
    Poor man's implementation of Chebyshev polynomial evaluation
    """

    assert degree >= 0

    if degree == 0:
        return 1
    elif degree == 1:
        return z
    else:
        return 2*z*chebyshev_poly_eval(degree-1,z) - chebyshev_poly_eval(degree-2,z)


def chebyshev_fit(x, y, N):
    
    M = len(x)  # panjang data
    assert M == len(y)

    z = transform_to_m1p1(x)
    
    # vektor b
    b = np.zeros(N+1)
    for i in range(N+1):
        s = 0.0
        for k in range(M):
            s = s + y[k]*chebyshev_poly_eval(i, z[k])
        b[i] = s

    A = np.zeros((N+1,N+1))
    for i in range(N+1):
        for j in range(N+1):
            s = 0.0
            for k in range(M):
                Ti = chebyshev_poly_eval(i, z[k])
                Tj = chebyshev_poly_eval(j, z[k])
                s = s + Ti*Tj
            A[i,j] = s
    print("A = ", A)
    print("b = ", b)
    c = np.linalg.solve(A,b)

    print("c = ", c)

    return c


def chebyshev_fit_eval(coefs, a, b, x):
    N = len(coefs) - 1
    z = (2*x - (a+b))/(b-a)
    s = 0.0
    for i in range(N+1):
        s = s + coefs[i]*chebyshev_poly_eval(i, z)
    return s


