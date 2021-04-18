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

def integ_romberg(f, a, b, es=1e-10, MAXIT=10):
    I = np.zeros( (MAXIT+2,MAXIT+2) )
    n = 1
    # We start from I[1,1], to follow the book's notation
    I[1,1] = integ_trapz_multiple(f, a, b, n)
    iterConv = 0
    for i in range(1,MAXIT+1):
        n = 2**i
        I[i+1,1] = integ_trapz_multiple(f, a, b, n)
        #
        for k in range(2,i+2):
            j = 2 + i - k
            I[j,k] = ( 4**(k-1)*I[j+1,k-1] - I[j,k-1] )/ (4**(k-1) - 1)
        #
        ea = abs( (I[1,i+1] - I[2,i])/I[1,i+1] )*100 # in percent
        if ea <= es:
            iterConv = i
            print("converged, iterConv = ", iterConv)
            break
        iterConv = i
    # to make sure that we are use variable that is defined outside the loop
    # we use iterConv instead of i
    return I[1,iterConv+1]
 




