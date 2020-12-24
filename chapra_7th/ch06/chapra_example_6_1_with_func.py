# Simple fixed-point iteration

# f(x) = exp(-x) - x
# The equation is transformed to
# x = exp(-x) = g(x)
#
# x_{i+1} = exp(-x_{i})

import numpy as np

def g(x):
    return np.exp(-x)


def root_fixed_point(g, x0, NiterMax=100, TOL=1e-10):
    # Initial guess
    x = x0
    for i in range(1,NiterMax+1):
        xnew = g(x)
        Δx = np.abs(xnew - x)
        print("%3d %18.10f %13.5e" % (i, xnew, Δx))
        # we are not using relative error here
        if Δx < TOL:
            break
        x = xnew
    return x

xroot = root_fixed_point(g, 0.0)
