import numpy as np

def optim_brent(f, xl, xu, TOL=1e-10, NiterMax=100, verbose=True):
    
    SMALL = np.finfo(np.float64).resolution

    ϕ = (1 + np.sqrt(5))/2
    ρ = 2.0 - ϕ
    u = xl + ρ*(xu - xl)
    v = u
    w = u
    x = u
    fu = f(u)
    fv = fu
    fw = fu
    fx = fu
    xm = 0.5*(xl + xu)
    d = 0.0
    e = 0.0

    iiter = 0

    while True:

        iiter = iiter + 1
        
        if verbose:
            print("\nBegin iter = ", iiter)
            print("x      = %18.10f" % x)
            print("xm     = %18.10f" % xm)
            print("fu     = %18.10f" % fu)
            print("x - xm = %18.10e" % abs(x - xm))

        if abs(x - xm) <= TOL:
            if verbose:
                print("Converged")
            break

        if iiter >= NiterMax:
            if verbose:
                print("WARNING: Maximum iterations reached")
            break

        para = abs(e) > TOL
        
        if para:
            # Try using parabolic interp
            r = (x - w)*(fx - fv)
            q = (x - v)*(fx - fw)
            p = (x - v)*q - (x - w)*r
            s = 2*(q - r)
            if s > 0:
                p = -q
            s = abs(s)
            #
            cond1 = abs(p) < abs(0.5*s*e)
            cond2 = p > s*(xl - x)
            cond3 = p < s*(xu - x)
            para = cond1 and cond2 and cond3
            # Parabolic interpolation step
            if para:
                if verbose:
                    print("Parabolic interpolation is used")
                e = d
                d = p/s

        if not para:
            if verbose:
                print("Using golden section")
            if x >= xm:
                e = xl - x
            else:
                e = xu - x
            d = ρ*e

        u = x + d
        fu = f(u)

        if fu <= fx:
            if u >= x:
                xl = x
            else:
                xu = x
            v = w; fv = fw
            w = x; fw = fx
            x = u; fx = fu
        else:
            #
            if u < x:
                xl = u
            else:
                xu = u
            #
            if (fu <= fw) or (abs(w - x) <= SMALL):
                v = w; fv = fw
                w = u; fw = fu
            elif (fu <= fv) or (abs(v - x) <= SMALL) or (abs(v - w) <= SMALL):
                v = u; fv = fu
        #
        xm = 0.5*(xl + xu)

    return xu, fu
