import numpy as np
from math import log, ceil

def bisection(f, x1, x2, TOL=1.0e-9, verbose=False, Niter=None):

    if verbose:
        print("Searching root with bisection method:")
        print("Interval = (%18.10f,%18.10f)" % (x1,x2))

    f1 = f(x1)
    if f1 == 0.0:
        return x1

    f2 = f(x2)
    if f2 == 0.0:
        return x2

    if np.sign(f1) == np.sign(f2):
        raise RuntimeError("Root is not bracketed")

    if Niter == None:
        Niter = int(ceil( log(abs(x2-x1)/TOL) )/ log(2.0) )
    
    if verbose:
        print("Bisection will iterate upto %d iterations" % (Niter))

    if verbose:
        print("Niter = ", Niter)

    # For the purpose of calculating relative error
    x3 = 0.0
    x3_old = 0.0

    for i in range(Niter):

        x3_old = x3
        x3 = 0.5*(x1 + x2)
        f3 = f(x3)

        if verbose:
            print("x3 = %18.10f" % (x3))

        if (abs(f3) > abs(f1)) and (abs(f3) > abs(f2)):
            return None

        if f3 == 0:
            return x3

        if np.sign(f2) != np.sign(f3):
            x1 = x3
            f1 = f3
        else:
            x2 = x3
            f2 = f3

    # return the result
    return x3, abs(x3 - x3_old)

