from math import sin, sqrt

def f(x):
    return sin(sqrt(x)) - x

# x = sin(sqrt(x))
def g(x):
    return sin(sqrt(x))

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

xroot = root_fixed_point(g, 0.5)
print("At root f(x) = ", f(xroot))
