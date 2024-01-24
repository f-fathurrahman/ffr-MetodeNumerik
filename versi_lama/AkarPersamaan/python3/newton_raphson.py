def newton_raphson(f, df, x, TOL=1.0e-10, verbose=False, NiterMax=100, multiplicity=1):

    if verbose:
        print("")
        print("Searching root with Newton-Raphson method:")
        print("Starting guess = %18.10f" % x)
        print("TOL = %e" % TOL)
        print("")

    assert TOL >= 0.0

    fx = f(x)
    if abs(fx) <= TOL:
        return x, 0.0

    x_old = 0.0

    for i in range(1,NiterMax+1):
        x_old = x

        dfx = df(x)
        if abs(dfx) <= TOL:
            raise RuntimeError("dfx is too small is newton_raphson")

        xroot = x - f(x)/dfx*multiplicity

        fxroot = f(xroot)

        if verbose:
            print("newton_raphson: %5d %18.10f %15.5e" % (i, xroot, abs(fxroot)))

        if abs(fxroot) <= TOL:
            print("Convergence is achived in %d iterations" % i)
            return xroot, abs(xroot-x_old)

        x = xroot

    return None, None
