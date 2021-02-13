from math import pi
from root_secant import *
from root_brent import *

R = 3.0  # m
V = 30.0 # m^3

def f(h):
    return V - pi*h**2*(3*R - h)/3

h = 3.0
#xroot = root_secant(f, h, h + 1.0)
xroot = root_brent(f, 0.0, h)
print(f(xroot))