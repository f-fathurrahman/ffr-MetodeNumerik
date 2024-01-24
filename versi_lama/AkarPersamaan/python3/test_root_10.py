import numpy as np
from fixed_point import *

def f(x):
    return np.exp(2*x) - 16*x**2

def g(x):
    return np.sqrt( np.exp(2*x)/16 )

x0 = 0.0
xroot, _ = fixed_point(g, x0, verbose=True)

print(xroot)
print(f(xroot))
print(f(-xroot))
