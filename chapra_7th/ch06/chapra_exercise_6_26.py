from math import sqrt
from root_fixed_point import *

scale = 1e-5
W = 1e6*scale
k = 0.25
V = 1e6*scale
Q = 1e5*scale

def f(c):
    return W - Q*c - k*V*sqrt(c)

def g1(c):
    return (W - Q*c/(k*V))**2

def g2(c):
    return (W - k*V*sqrt(c))/Q

#xr = root_fixed_point(g1, 3.0, NiterMax=5)
xr = root_fixed_point(g2, 3.0, NiterMax=100)
print("At root: ", f(xr))
#xr = root_fixed_point(g1, 3.0, NiterMax=10)