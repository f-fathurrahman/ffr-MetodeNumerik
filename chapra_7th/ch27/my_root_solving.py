from math import ceil, log10

def root_bisection(f, x1, x2, TOL=1.0e-9, NiterMax=None ):

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

    print(13*" "+"Iter      Estimated          f(x)")
    print(13*" "+"----      ---------          ----")
    print("")

    for i in range(1,NiterMax+1):

        x3_old = x3
        x3 = 0.5*(x1 + x2)
        f3 = f(x3)

        print("bisection: %5d %18.10f %15.5e" % (i, x3, abs(f3)))

        if abs(f3) <= TOL:
            print("")
            print("bisection is converged in %d iterations" % i)
            # return the result
            return x3

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

    print("No root is found")
    return None

def root_regula_falsi(f, xl, xu, NiterMax=100, TOL=1e-10):
    fxl = f(xl)
    fxu = f(xu)
    assert xl < xu
    assert fxl*fxu < 0
    #
    for i in range(1,NiterMax+1):
        xr = xu - fxu*(xl - xu)/(fxl - fxu)
        fxr = f(xr)
        print("regula_falsi: %4d %18.10f %18.10e" % (i, xr, fxr))
        if abs(fxr) < TOL:
            break # convergence achieved
        if fxl*fxr < 0:
            xu = xr
            fxu = fxr
        else:
            xl = xr
            fxl = fxr
    return xr
