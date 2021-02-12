from math import exp
from root_secant import *
from root_newton_raphson import *

def f(x):
    return x**3 - 2*x**2 - 4*x + 8

def df(x):
    return 3*x**2 - 4*x - 4 

def d2f(x):
    return 6*x - 4

# Guess
x = 1.2

#xroot = root_newton_raphson(f, df, x)
xroot = root_newton_raphson_mod(f, df, d2f, x)
#xroot = root_newton_raphson(f, df, x, m=2)
#xroot = root_secant(f, x, x + 1.0)
#xroot = root_secant_mod(f, x)
print(f(xroot))