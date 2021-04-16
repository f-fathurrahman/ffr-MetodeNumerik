import numpy as np

def integ_trapz_multiple( f, a, b, N ):
    x0 = a
    xN = b
    h = (b - a)/N
    ss = 0.0
    for i in range(1,N):
        xi = x0 + i*h
        ss = ss + f(xi)
    I = h/2 * ( f(x0) + 2*ss + f(xN) )
    return I


def integ_trapz( f, a, b ):
    I = (b - a) * (f(a) + f(b)) / 2
    return I

def integ_trapz_multiple_v2(f, a, b, N):
    s = 0.0
    Δ = (b-a)/N
    for i in range(N):
        aa = a + i*Δ
        bb = a + (i+1)*Δ
        s = s + integ_trapz(f, aa, bb)
    return s


# N is number of intervals
# N+1 is number of points (total function evaluations)
def integ_simpson13_multiple( f, a, b, N ):
    
    assert N >= 2

    x0 = a
    xN = b
    h = (b - a)/N
    x = np.linspace(a, b, N+1)
    
    ss_odd = 0.0
    for i in range(1,N,2):
        ss_odd = ss_odd + f(x[i])
    
    ss_even = 0.0
    for i in range(2,N-1,2):
        ss_even = ss_even + f(x[i])

    I = (b - a)/(3*N) * ( f(x0) + 4*ss_odd + 2*ss_even + f(xN) )
    return I

def integ_simpson13( f, a, b ):
    #
    h = (b - a)/2
    x0 = a
    x1 = a + h
    x2 = b
    #
    I = h/3 * ( f(x0) + 4*f(x1) + f(x2) )
    return I

def integ_simpson13_multiple_v2(f, a, b, N):
    s = 0.0
    Δ = (b-a)/N
    for i in range(N):
        aa = a + i*Δ
        bb = a + (i+1)*Δ
        s = s + integ_simpson13(f, aa, bb)
    return s

def integ_simpson38( f, a, b ):
    #
    h = (b - a)/3
    x0 = a
    x1 = a + h
    x2 = a + 2*h
    x3 = b
    #
    I = 3*h/8 * ( f(x0) + 3*f(x1) + 3*f(x2) + f(x3) )
    return I

def integ_simpson38_multiple(f, a, b, N):
    s = 0.0
    Δ = (b-a)/N
    for i in range(N):
        aa = a + i*Δ
        bb = a + (i+1)*Δ
        s = s + integ_simpson38(f, aa, bb)
    return s

def integ_boole(f, a, b):
    h = (b-a)/4
    x0 = a
    x1 = a + h
    x2 = a + 2*h
    x3 = a + 3*h
    x4 = b

    s = 7*f(x0) + 32*f(x1) + 12*f(x2) + 32*f(x3) + 7*f(x4)
    return (b-a)*s/90

def integ_boole_multiple(f, a, b, N):
    s = 0.0
    Δ = (b-a)/N
    for i in range(N):
        aa = a + i*Δ
        bb = a + (i+1)*Δ
        s = s + integ_boole(f, aa, bb)
    return s