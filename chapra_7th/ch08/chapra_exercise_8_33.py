from math import sqrt
from my_root_finders import root_bisection

Z = 75.0
R = 225.0
C = 0.6e-6
L = 0.5

def f(ω):
    return 1/Z - sqrt( 1/R**2 + (ω*C - 1/(ω*L))**2 )

xroot = root_bisection(f, 1.0, 1000.0)
print("xroot = ", xroot)
print("Check function value at root: ", f(xroot))
