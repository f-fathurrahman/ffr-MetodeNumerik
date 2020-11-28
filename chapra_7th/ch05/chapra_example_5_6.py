import numpy as np

def f(x):
    return x**10 - 1

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
            ε_a = abs((xr - xr_old)/xr)*100
            print("%10.1f" % ε_a)
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
            ε_a = abs((xr - xr_old)/xr)*100
            print("%10.1f" % ε_a)
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


# Initial guess
xl = 0.0
xu = 1.3
print("Using bisection:")
root_bisection(f, xl, xu, NiterMax=5)
#root_regula_falsi(f, xl, xu, NiterMax=200)

print("\nUsing regula falsi:")
root_regula_falsi(f, xl, xu, NiterMax=5)
#root_regula_falsi(f, xl, xu, NiterMax=200)