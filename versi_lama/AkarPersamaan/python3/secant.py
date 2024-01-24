def secant(f, x, TOL=1.0e-10, verbose=False, NiterMax=100, multiplicity=1):

    if verbose:
        print("")
        print("Searching root with secant method:")
        print("Starting guess = %18.10f" % x)
        print("TOL = %e" % TOL)
        print("")

    assert TOL >= 0.0

    fx = f(x)
    if abs(fx) <= TOL:
        return x, 0.0

    x_old = 0.0

    for i in range(1,NiterMax+1):

        fx = f(x)

        if i == 1:
            SMALL = 1e-10
            fx_old = f(x+SMALL)
            dfx = abs(fx - fx_old)/SMALL
        else:
            dfx = (fx - fx_old)/(x - x_old)

        if abs(dfx) <= TOL:
            raise RuntimeError("dfx is too small in secant")

        x_old = x
        fx_old = fx

        xroot = x - f(x)/dfx*multiplicity

        fxroot = f(xroot)

        if verbose:
            print("secant: %5d %18.10f %15.5e" % (i, xroot, abs(fxroot)))

        if abs(fxroot) <= TOL:
            print("Convergence is achived in %d iterations" % i)
            return xroot, abs(xroot-x_old)

        x = xroot

    return None, None
