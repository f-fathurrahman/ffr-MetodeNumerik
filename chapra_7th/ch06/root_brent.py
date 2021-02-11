import numpy as np

def root_brent(f, xl, xu, TOL=1.0e-9, verbose=False, NiterMax=100):

    if verbose:
        print("")
        print("Searching root with Brent's method:")
        print("Interval = (%18.10f,%18.10f)" % (xl,xu))
        print("TOL = %e" % TOL)
        print("")

    a = xl
    b = xu
    fa = f(a)
    fb = f(b)

    c = a
    fc = fa
    d = b - c
    e = d

    b_old = 0.0

    for i in range(1,NiterMax+1):

        b_old = b

        if abs(fb) == 0:
            return b

        if fa*fb > 0:
            a = c
            fa = fc
            d = b - c
            e = d

        if abs(fa) < abs(fb):
            c = b
            b = a
            a = c

            fc = fb
            fb = fa
            fa = fc

        m = 0.5*(a - b)
        if abs(m) <= TOL or abs(fb) <= TOL:
            break

        if abs(e) >= TOL and abs(fc) > abs(fb):
            s = fb/fc
            if a == c:
                p = 2*m*s
                q = 1 - s
            else:
                q = fc/fa
                r = fb/fa
                p = s*( 2*m*q*(q-r) - (b-c)*(r-1) )
                q = (q-1)*(r-1)*(s-1)

            if p > 0.0:
                q = -q
            else:
                p = -p

            if 2*p < (3*m*q - abs(TOL*q)) and p < abs(0.5*e*q):
                e = d
                d = p/q
            else:
                d = m
                e = m
        else:
            d = m
            e = m
            
        c = b
        fc = fb

        if abs(d) > TOL:
            b = b + d
        else:
            b = b - np.sign(b-a)*TOL
    
        fb = f(b)

        if verbose:
            print("brent: %5d %18.10f %15.5e" % (i, b, abs(fb)))

    return b, abs(b - b_old)

