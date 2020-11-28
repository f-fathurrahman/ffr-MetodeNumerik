import numpy as np

def f(x):
    return x**3.5 - 80

def root_bisection(f, xl, xu, NiterMax=100, eps_a=1e-10):
    fxl = f(xl)
    fxu = f(xu)
    assert xl < xu
    assert fxl*fxu < 0
    ε_a = 1e10 # set to big number
    SMALL = np.finfo(float).eps
    #
    for i in range(1,NiterMax+1):
        xr = 0.5*(xl + xu)
        fxr = f(xr)
        print("%4d %18.10f %18.10e" % (i, xr, fxr), end="")
        if i > 1 and abs(xr) > SMALL:
            ε_a = abs((xr - xr_old)/xr)
            print("%10.1f%%" % (ε_a*100))
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

# Part c
xl = 2.0
xu = 5.0
print("Using bisection: ")
root_bisection(f, xl, xu, eps_a=2.5/100)

print("Using regula falsi: ")
root_regula_falsi(f, xl, xu, eps_a=2.5/100)

