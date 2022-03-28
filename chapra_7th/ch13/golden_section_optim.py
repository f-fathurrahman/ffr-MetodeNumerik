import numpy as np

# search for maximum value

def golden_section_optim(f, xlow, xhigh, NiterMax=100, TOL=1e-10):

    SMALL = np.finfo(np.float64).resolution

    γ = (np.sqrt(5) - 1)/2
    xl = xlow
    xu = xhigh
    iiter = 1
    d = γ*(xu - xl)
    x1 = xl + d
    x2 = xu - d
    f1 = f(x1)
    f2 = f(x2)
    
    fx_old = np.nan # use this for evaluating convergence?

    if f1 > f2:
        xopt = x1
        fx = f1
    else:
        xopt = x2
        fx = f2

    for iiter in range(NiterMax):
        d = γ*d
        xint = xu - xl
        if f1 > f2:
            xl = x2
            x2 = x1
            x1 = xl + d
            f2 = f1
            f1 = f(x1)
        else:
            xu = x1
            x1 = x2
            x2 = xu - d
            f1 = f2
            f2 = f(x2)
        # Update iiter?
        if f1 > f2:
            xopt = x1
            fx = f1
        else:
            xopt = x2
            fx = f2
        #
        if abs(xopt) > SMALL:
            ea = (1 - γ)*abs(xint/xopt)

        print("%3d %18.10f %18.10f %15.5e %15.5e" % (iiter, xopt, fx, ea, d))

        if ea <= TOL:
            print("Converged")
            break

    return xopt, fx

