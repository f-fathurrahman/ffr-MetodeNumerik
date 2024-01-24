
"""
Search for root of the equation: x = g(x) using fixed-point iteration.
"""
def fixed_point(g, x, TOL=1.0e-10, verbose=False, NiterMax=100, multiplicity=1):

    if verbose:
        print("")
        print("Searching root with fixed-point method:")
        print("Starting guess = %18.10f" % x)
        print("TOL = %e" % TOL)
        print("")

    assert TOL >= 0.0

    gx = g(x)
    if abs(x - gx) <= TOL:
        return x, 0.0

    x_old = 0.0

    for i in range(1,NiterMax+1):
        
        x_old = x
        x = gx

        gx = g(x)

        if verbose:
            print("fixed_point %5d %18.10f %15.5e" % (i, x, abs(x-gx)))

        if abs(x-gx) <= TOL:
            print("Convergence is achived in %d iterations" % i)
            return x, abs(x-x_old)

    return None, None
