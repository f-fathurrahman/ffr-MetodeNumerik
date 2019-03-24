import numpy as np
from math import sin, cos, pi

# x must have homogeneous interval
def fourier_fit(x, y, m=1):
    
    N = len(x)

    assert N == len(y)

    assert N > (2*m + 1)

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

    A = 2.0*A/N
    B = 2.0*B/N

    return A0, A, B


def fourier_fit_eval(A0, A, B, x, xo):

    m = len(A)
    N = len(x)

    xmin = np.min(x)
    xmax = np.max(x)
    s = (x - xmin)*(N-1)*2*pi/(N*(xmax-xmin))
    so = (xo - xmin)*(N-1)*2*pi/(N*(xmax-xmin))
    y = A0
    for i in range(1,m+1):
        y = y + A[i-1]*cos(i*so) + B[i-1]*cos(i*so)
    return y


def sinusoid_fit(x, y, equispaced=True, T=0.0):
    N = len(x)

    assert N == len(y)

    if equispaced:
        deltax = x[1] - x[0]
        T = (N-1)*deltax
        print("deltax = ", deltax)
    else:
        if T <= 0:
            raise RuntimeError("Perlu T > 0")
    
    omega = 2*pi/T

    if equispaced:

        A0 = np.sum(y)/N

        A1 = 0.0
        B1 = 0.0
        for i in range(N):
            A1 = A1 + y[i]*cos(omega*x[i])
            B1 = B1 + y[i]*sin(omega*x[i])
        A1 = 2.0*A1/N
        B1 = 2.0*B1/N

        return A0, A1, B1, omega

    else:

        NormalMatrix = np.zeros((3,3))
        RHSVector = np.zeros(3)

        NormalMatrix[0,0] = N
        NormalMatrix[1,1] = np.sum( np.cos(omega*x)**2 )
        NormalMatrix[2,2] = np.sum( np.sin(omega*x)**2 )

        NormalMatrix[0,1] = np.sum( np.cos(omega*x) )
        NormalMatrix[1,0] = NormalMatrix[0,1]

        NormalMatrix[0,2] = np.sum( np.sin(omega*x) )
        NormalMatrix[2,0] = NormalMatrix[0,2]

        NormalMatrix[2,1] = np.sum( np.cos(omega*x)*np.sin(omega*x) )
        NormalMatrix[1,2] = NormalMatrix[2,1]

        RHSVector[0] = np.sum(y)
        RHSVector[1] = np.sum(y*np.cos(omega*x))
        RHSVector[2] = np.sum(y*np.sin(omega*x))

        sol = np.linalg.solve(NormalMatrix, RHSVector)

        return sol[0], sol[1], sol[2], omega



def sinusoid_fit_eval(A0, A1, B1, omega, x):
    return A0 + A1*cos(omega*x) + B1*sin(omega*x)

