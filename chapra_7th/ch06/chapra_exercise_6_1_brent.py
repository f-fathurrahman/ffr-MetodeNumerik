from math import sin, sqrt
from root_brent import *

def f(x):
    return sin(sqrt(x)) - x

def root_fixed_point(g, x0, NiterMax=100, TOL=1e-10):
    # Initial guess
    x = x0
    for i in range(1,NiterMax+1):
        xnew = g(x)
        Δx = abs(xnew - x)
        print("%3d %18.10f %13.5e" % (i, xnew, Δx))
        # we are not using relative error here
        if Δx < TOL:
            x = xnew
            break
        x = xnew
    return x

xroot, _ = root_brent(f, 0.5, 1.0, verbose=True)
print("At root f(x) = ", f(xroot))
