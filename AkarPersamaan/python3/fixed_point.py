def fixed_point(f, x, TOL=1.0e-10, verbose=False, NiterMax=100, multiplicity=1):

    if verbose:
        print("")
        print("Searching root with fixed-point method:")
        print("Starting guess = %18.10f" % x)
        print("TOL = %e" % TOL)
        print("")

    assert TOL >= 0.0

    fx = f(x)
    if abs(x - fx) <= TOL:
        return x, 0.0

    x_old = 0.0

    for i in range(1,NiterMax+1):
        
        x_old = x
        x = fx

        fx = f(x)

        if verbose:
            print("fixed_point %5d %18.10f %15.5e" % (i, x, abs(x-fx)))

        if abs(x-fx) <= TOL:
            print("Convergence is achived in %d iterations" % i)
            return x, abs(x-x_old)

    return None, None
