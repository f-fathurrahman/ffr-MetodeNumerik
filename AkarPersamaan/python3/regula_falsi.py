def regula_falsi(f, x1, x2, TOL=1.0e-9, verbose=False, NiterMax=100):

    if verbose:
        print("")
        print("Searching root with regula falsi method:")
        print("Interval = (%18.10f,%18.10f)" % (x1,x2))
        print("TOL = %e" % TOL)        
        print("")

    assert TOL >= 0.0

    f1 = f(x1)
    if abs(f1) <= TOL:
        return x1, 0.0

    f2 = f(x2)
    if abs(f2) <= TOL:
        return x2, 0.0

    if f1*f2 > 0.0:
        raise RuntimeError("Root is not bracketed")    

    # For the purpose of calculating relative error
    x3 = 0.0
    x3_old = 0.0

    for i in range(1,NiterMax+1):

        x3_old = x3

        x3 = (x1*f2 - x2*f1)/(f2 - f1)
        f3 = f(x3)
        
        if verbose:
            print("regula_falsi: %5d %18.10f %18.10f" % (i, x3, f3))

        if abs(f3) <= TOL:
            if verbose:
                print("")
                print("regula_falsi: Convergence is achieved in %d iterations" % (i))
            break

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

    return x3, abs(x3 - x3_old)
