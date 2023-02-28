import numpy as np

# This will search for minimum value
# If you want to use this for searching maximum value please
# define f = -f_actual, where f_actual is the function you want
# to maximize.
# The sign keyword argument is not implemented here.

def optim_golden_ratio(f, xlow, xhigh, NiterMax=100, TOL=1e-10, verbose=True):

    SMALL = np.finfo(np.float64).resolution

    ϕ = (np.sqrt(5) - 1)/2
    xl = xlow
    xu = xhigh
    iiter = 1
    d = ϕ*(xu - xl)
    x1 = xl + d
    x2 = xu - d
    f1 = f(x1)
    f2 = f(x2)

    if verbose:
        print("Iteration: ", 0)
        print()
        print("xl = %18.10f xu = %18.10f" % (xl, xu))
        print("-"*47)
        print("x2 = %18.10f x1 = %18.10f" % (x2, x1))
        print("f2 = %18.10f f1 = %18.10f" % (f2, f1))

    if f1 < f2:
        xopt = x1
        fxopt = f1
        if verbose:
            print("-"*47)
            print("f1 is the current optimum value")
            print("(xopt,fxopt) = (%18.10f,%18.10f)" % (xopt, fxopt))
    else:
        xopt = x2
        fxopt = f2
        if verbose:
            print("-"*47)
            print("f2 is the current optimum value")
            print("(xopt,fxopt) = (%18.10f,%18.10f)" % (xopt, fxopt))

    for iiter in range(1,NiterMax+1):

        d = ϕ*d
        xint = xu - xl

        if f1 < f2:
            xl = x2
            x2 = x1
            x1 = xl + d
            f2 = f1
            f1 = f(x1)
            if verbose:
                print("")
                print("For next iteration: ")
                print("- Replacing xl with x2")
                print("- New point: x1 = %18.10f" % x1)
                print("- xu is not changed")
                print("- Replacing x2 with x1")
        else:
            xu = x1
            x1 = x2
            x2 = xu - d
            f1 = f2
            f2 = f(x2)
            if verbose:
                print("")
                print("For next iteration: ")
                print("- xl is not changed")
                print("- x1 is replaced by x2")
                print("- New point: x2 = %18.10f" % x2)
                print("- xu is replaced by x1")

        if verbose:
            print("")
            print("Iteration: ", iiter)
            print("xl = %18.10f xu = %18.10f" % (xl, xu))
            print("-"*47)
            print("x2 = %18.10f x1 = %18.10f" % (x2, x1))
            print("f2 = %18.10f f1 = %18.10f" % (f2, f1))

        #
        if f1 < f2:
            xopt = x1
            fxopt = f1
            if verbose:
                print("-"*47)
                print("f1 is the current optimum value")
                print("(xopt,fxopt) = (%18.10f,%18.10f)" % (xopt, fxopt))
        else:
            xopt = x2
            fxopt = f2
            if verbose:
                print("-"*47)
                print("f2 is the current optimum value")
                print("(xopt,fxopt) = (%18.10f,%18.10f)" % (xopt, fxopt))

        if abs(xopt) > SMALL:
            ea = (1 - ϕ)*abs(xint/xopt)
        else:
            # The above might fail if xopt is very close to zero
            # We set xint as the convergence criteria
            ea = xint

        if verbose:
            print("Interval length = %18.10e" % xint)
            print("ea              = %18.10e" % ea)

        if ea <= TOL:
            if verbose:
                print("Converged")
            break

    return xopt, fxopt

