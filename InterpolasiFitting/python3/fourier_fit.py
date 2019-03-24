import numpy as np
from math import sin, cos, pi

# x must have homogeneous interval
def fourier_fit(x, y, m=1):
    
    N = len(x)

    assert N == len(y)

    assert N > 2*m + 1

    s = np.zeros(N)
    for i in range(N):
        s[i] = i*2*pi/N

    A0 = 1/N*np.sum(y)

    A = np.zeros(m)  # coef for cos
    B = np.zeros(m)  # coef for sin

    for i in range(1,m+1):
        for j in range(N):
            A[i-1] = A[i-1] + y[j]*cos(i*s[j])
            B[i-1] = B[i-1] + y[j]*sin(i*s[j])

    A = 2*A/N
    B = 2*B/N

    return A0, A, B


def fourier_fit_eval(A0, A, B, x, xo):

    m = len(A)
    N = len(x)

    xmin = np.min(x)
    xmax = np.max(x)
    s = (x - xmin)*(N-1)*2*pi/(N*(xmax-xmin))
    so = (xo - xmin)*(N-1)*2*pi/(N*(xmax-xmin))
    print("so = ", so)
    y = A0
    for i in range(1,m+1):
        y = y + A[i-1]*cos(i*so) + B[i-1]*cos(i*so)
    return y


