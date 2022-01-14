from math import sqrt, pow
import cmath
import numpy as np

def root_bisection(f, xl, xu, NiterMax=100, TOL=1e-10, verbose=False):
    fxl = f(xl)
    fxu = f(xu)
    assert xl < xu
    assert fxl*fxu < 0
    #
    for i in range(1,NiterMax+1):
        xr = 0.5*(xl + xu)
        fxr = f(xr)
        if verbose:
            print("bisection: %4d %18.10f %18.10e" % (i, xr, fxr))
        if abs(fxr) < TOL:
            print("Convergence achived in %d iterations" % i)
            break # convergence achieved
        if fxl*fxr < 0:
            xu = xr
            fxu = fxr
        else:
            xl = xr
            fxl = fxr
    return xr


def root_brent(f, xl, xu, TOL=1.0e-10, NiterMax=100, verbose=True):
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

    return b



def root_fixed_point(g, x0, NiterMax=100, TOL=1e-10):
    # Initial guess
    x = x0
    for i in range(1,NiterMax+1):
        xnew = g(x)
        Δx = abs(xnew - x)
        print("fixed_point: %3d %18.10f %13.5e" % (i, xnew, Δx))
        # we are not using relative error here
        if Δx < TOL:
            x = xnew
            break
        x = xnew
    return x


def root_muller(f, xr, h=0.5, TOL=1.0e-10, NiterMax=100):

    res = 0.0;
    a = xr - h
    b = xr
    c = xr + h

    for i in range(1,NiterMax+1): 
      
        # Calculating various constants  
        # required to calculate x3 
        f1 = f(a)
        f2 = f(b)
        f3 = f(c)
        
        d1 = f1 - f3
        d2 = f2 - f3
        
        h1 = a - c
        h2 = b - c
        
        a0 = f3
        a1 = (((d2 * h1**2) - (d1 * h2**2)) / ((h1 * h2) * (h1 - h2))) 
        a2 = (((d1 * h2) - (d2 * h1)) / ( h1*h2*(h1 - h2) ))
        discr = a1 * a1 - 4 * a0 * a2
        is_real = discr >= 0.0
        if is_real:
            x = ((-2 * a0) / (a1 + sqrt(discr) ))
            y = ((-2 * a0) / (a1 - sqrt(discr) ))
        else:
            x = ((-2 * a0) / (a1 + cmath.sqrt(discr) ))
            y = ((-2 * a0) / (a1 - cmath.sqrt(discr) ))

        # Taking the root which is closer to x2
        if is_real:
            if x >= y: 
                res = x + c
            else: 
                res = y + c 
        else:
            # FIXME: need to check this!!!
            if abs(x-c) >= abs(y-c): 
                res = x + c
            else: 
                res = y + c 

  
        print("muller: %3d %18.10f %18.10e" % (i, res, abs(f(res))))
        if abs(f(res)) <= TOL:
            break
        
        res_old = res
        a = b; 
        b = c; 
        c = res;

    return res

# m is root multiplicity
def root_newton_raphson(f, df, x0, NiterMax=100, TOL=1e-10, m=1):
    x = x0
    for i in range(1,NiterMax+1):
        xnew = x  - m*f(x)/df(x) # Newton-Raphson formula
        fxnew = f(xnew)
        print("newton_raphson: %3d %18.10f %18.10e" % (i, xnew, fxnew))
        if abs(f(xnew)) < TOL:
            x = xnew
            break
        x = xnew
    return x

# Modified Newton-Raphson
# Need 2nd derivative info of f
def root_newton_raphson_mod(f, df, d2f, x0, NiterMax=100, TOL=1e-10):
    x = x0
    for i in range(1,NiterMax+1):
        fx = f(x)
        dfx = df(x)
        d2fx = d2f(x)
        # Using equation 6.16
        xnew = x  - fx*dfx/( dfx**2 - fx*d2fx )
        fxnew = f(xnew)
        print("newton_raphson_mod: %3d %18.10f %18.10e" % (i, xnew, fxnew))
        if abs(f(xnew)) < TOL:
            x = xnew
            break
        x = xnew
    return x

def root_regula_falsi_mod(f, xl, xu, NiterMax=100, TOL=1e-10):
    fxl = f(xl)
    fxu = f(xu)
    assert xl < xu
    assert fxl*fxu < 0
    #
    il = 0
    iu = 0
    for i in range(1,NiterMax+1):
        xr = xu - fxu*(xl - xu)/(fxl - fxu)
        fxr = f(xr)
        print("regula_falsi_mod: %4d %18.10f %18.10e" % (i, xr, fxr))
        if abs(fxr) < TOL:
            break # convergence achieved
        if fxl*fxr < 0:
            xu = xr
            fxu = fxu
            iu = 0
            il = il + 1
            if il >= 2:
                fxl = 0.5*fxl
        else:
            xl = xr
            fxl = fxr
            il = 0
            iu = iu + 1
            if iu >= 2:
                fxu = 0.5*fxu
    return xr

def root_regula_falsi(f, xl, xu, NiterMax=100, TOL=1e-10):
    fxl = f(xl)
    fxu = f(xu)
    assert xl < xu
    assert fxl*fxu < 0
    #
    for i in range(1,NiterMax+1):
        xr = xu - fxu*(xl - xu)/(fxl - fxu)
        fxr = f(xr)
        print("regula_falsi: %4d %18.10f %18.10e" % (i, xr, fxr))
        if abs(fxr) < TOL:
            break # convergence achieved
        if fxl*fxr < 0:
            xu = xr
            fxu = fxr
        else:
            xl = xr
            fxl = fxr
    return xr



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


def root_secant(f, x0, x1, NiterMax=100, TOL=1e-10):
    for i in range(1,NiterMax+1):
        # approximation of derivative of f(x)
        dfx = (f(x0) - f(x1))/(x0 - x1)
        if abs(dfx) < 1e-12:
            print("WARNING: small dfx = ", dfx)
        xnew = x1 - f(x1)/dfx
        #
        fxnew = f(xnew)
        print("secant: %3d %18.10f %18.10e" % (i, xnew, fxnew))
        if abs(f(xnew)) < TOL:
            x1 = xnew
            break
        x0 = x1
        x1 = xnew
    return x1

def root_secant_mod(f, x1, δ=0.01, NiterMax=100, TOL=1e-10):
    for i in range(1,NiterMax+1):
        # approximation of derivative of f(x)
        dfx = (f(x1+δ) - f(x1))/δ
        if abs(dfx) < 1e-12:
            print("WARNING: small dfx = ", dfx)
        xnew = x1 - f(x1)/dfx
        #
        fxnew = f(xnew)
        print("%3d %18.10f %18.10e" % (i, xnew, fxnew))
        if abs(f(xnew)) < TOL:
            x1 = xnew
            break
        x1 = xnew
    return x1
