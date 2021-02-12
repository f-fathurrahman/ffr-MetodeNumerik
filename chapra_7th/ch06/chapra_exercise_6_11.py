from math import exp
from root_secant import *
from root_newton_raphson import *

def f(x):
    return exp(-0.5*x)*(4 - x) - 2

def df(x):
    return -0.5*(4 - x)*exp(-0.5*x) - exp(-0.5*x)

# Guess
#x = 2.0
x = 6.0 # div by zero
#x = 8.0

xroot = root_newton_raphson(f, df, x)
#xroot = root_secant(f, x, x + 1.0)
#xroot = root_secant_mod(f, x)
print(f(xroot))