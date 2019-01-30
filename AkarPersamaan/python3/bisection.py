from math import log10, ceil

def bisection(f, x1, x2, TOL=1.0e-9, verbose=False, NiterMax=None ):

    if verbose:
        print("")
        print("Searching root with bisection method:")
        print("Interval = (%18.10f,%18.10f)" % (x1,x2))
        print("TOL = %e" % TOL)
        print("")

    f1 = f(x1)
    if abs(f1) <= TOL:
        return x1, 0.0

    f2 = f(x2)
    if abs(f2) <= TOL:
        return x2, 0.0

    if f1*f2 > 0.0:
        raise RuntimeError("Root is not bracketed")

    # No NiterMax is provided
    # We calculate the default value here.
    if NiterMax == None:
        NiterMax = int(ceil( log10(abs(x2-x1)/TOL) )/ log10(2.0) ) + 10
        # extra 10 iterations

    # For the purpose of calculating relative error
    x3 = 0.0
    x3_old = 0.0

    if verbose:
        print(13*" "+"Iter      Estimated          f(x)")
        print(13*" "+"----      ---------          ----")
        print("")


    for i in range(1,NiterMax+1):

        x3_old = x3
        x3 = 0.5*(x1 + x2)
        f3 = f(x3)

        if verbose:
            print("bisection: %5d %18.10f %15.5e" % (i, x3, abs(f3)))

        if abs(f3) <= TOL:
            if verbose:
                print("")
                print("bisection is converged in %d iterations" % i)
            # return the result
            return x3, abs(x3 - x3_old)


        if f2*f3 < 0.0:
            # sign of f2 and f3 is different
            # root is in [x2,x3]
            # change the interval bound of x1 to x3
            x1 = x3
            f1 = f3
        else:
            # sign of f1 and f3 is different
            # root is in [x1,x3]
            # change the interval bound of x2 to x3
            x2 = x3
            f2 = f3


    # No root is found after NiterMax iterations
    if verbose:
        print("No root is found")
    return None, None

