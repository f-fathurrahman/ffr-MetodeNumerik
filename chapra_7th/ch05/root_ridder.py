from math import sqrt

def root_ridder(f, x1, x2, NiterMax=100, TOL=1.0e-10):

    assert TOL >= 0.0

    f1 = f(x1)
    if abs(f1) <= TOL:
        return x1

    f2 = f(x2)
    if abs(f2) <= TOL:
        return x2

    if f1*f2 > 0.0:
        raise RuntimeError("Root is not bracketed")

    # For the purpose of calculating relative error
    x3 = 0.0

    for i in range(1,NiterMax+1):

        c = 0.5*(x1 + x2)
        fc = f(c)

        s = sqrt(fc**2 - f1*f2)
        if s == 0.0:
            raise RuntimeError("s is zero in ridder")

        dx = (c - x1)*fc/s
        if (f1 - f2) < 0:
            dx = -dx

        # new approximation of root
        x3 = c + dx
        f3 = f(x3)
        
        print("ridder: %5d %18.10f %15.5e" % (i, x3, abs(f3)))

        if abs(f3) <= TOL:
            return x3

        # Rebracket
        if fc*f3 > 0.0:
            if f1*f3 < 0.0:
                x2 = x3
                f2 = f3
            else:
                x1 = x3
                f1 = f3
        else:
            x1 = c
            x2 = x3
            f1 = fc
            f2 = f3

    # No root is found after NiterMax iterations
    print("No root is found returning last value")
    return x3
