import numpy as np

def f(x):
    return np.log(x)

def root_regula_falsi(f, xl, xu, NiterMax=100, eps_a=1e-10):
    fxl = f(xl)
    fxu = f(xu)
    assert xl < xu
    assert fxl*fxu < 0
    ε_a = 1e10 # set to big number
    SMALL = np.finfo(float).eps
    #
    for i in range(1,NiterMax+1):
        xr = xu - fxu*(xl - xu)/(fxl - fxu)
        fxr = f(xr)
        print("%4d %18.10f %18.10e" % (i, xr, fxr), end="")
        if i > 1 and abs(xr) > SMALL:
            ε_a = abs((xr - xr_old)/xr)
            print("%10.1f%%" % (ε_a*100)) # in percent
        else:
            print("")
        if ε_a < eps_a:
            break # convergence achieved
        xr_old = xr
        if fxl*fxr < 0:
            xu = xr
            fxu = fxr
        else:
            xl = xr
            fxl = fxr
    return xr


#def root_secant(f, x0, x1)
#def root_secant(f, x0, x1, NiterMax=100, TOL=1e-10):


x0 = 0.5
x1 = 5.0
for i in range(1,4):
    # approximation of derivative of f(x)
    dfx = (f(x0) - f(x1))/(x0 - x1)
    #
    xnew = x1 - f(x1)/dfx
    fxnew = f(xnew)
    print("%3d %18.10f %18.10e" % (i, xnew, fxnew))
    x0 = x1
    x1 = xnew

x0 = 0.5
x1 = 5.0
xr = root_regula_falsi(f, x0, x1)