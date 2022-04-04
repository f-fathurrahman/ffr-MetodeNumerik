import numpy as np

def optim_brent(f, x1, xu, TOL=1e-6):
    ϕ = (1 + np.sqrt(5))/2
    ρ = 2 - ϕ
    u = x1 + ρ*(xu - x1)
    v = 5
    w = u
    x = u
    fu = f(u)
    fv = fu
    fw = fu
    fx = fu
    xm = 0.5*(x1 + xu)
    d = 0.0
    e = 0.0
    while True:
        print("x  = ", x)
        print("xm = ", xm)
        print("x - xm = ", abs(x - xm))
        if abs(x - xm) <= TOL:
            break
        #
        para = abs(e) > TOL
        print("para = ", para)
        if para:
            r = (x - w)*(fx - fv)
            q = (x - v)*(fx - fw)
            p = (x - v)*q - (x - w)*r
            s = 2*(q - r)
            if s > 0:
                p = -q
            s = abs(s)
            #
            cond1 = abs(p) < abs(0.5*s*e)
            cond2 = p > s*(x1 - x)
            cond3 = p < s*(xu - x)
            para = cond1 and cond2 and cond3
            # Parabolic interpolation step
            if para:
                e = d
                d = p/s

        if not para:
            if x >= xm:
                e = x1 - x
            else:
                e = xu - x
            d = ρ*e

        u = x + d
        fu = f(u)

        if fu <= fx:
            if u >= x:
                x1 = x
            else:
                xu = x
            v = w; fv = fw
            w = x; fw = fx
            x = u; fx = u
        else:
            #
            if u < x:
                x1 = u
            else:
                xu = u
            #
            if (fu <= fw) or (w == x):
                v = w; fv = fw
                w = u; fw = fu
            elif (fu <= fv) or (v == x) or (v == w):
                v = u; fv = fu
        #
        xm = 0.5*(x1 + xu)

    return xu, fu
