from bisection import *

def f(x):
    return x**3 - 10*x**2 + 5.0

x1 = 0.0
x2 = 1.0

x = bisection(f, x1, x2, TOL=1e-4)
print("x = {:6.4f}".format(x))