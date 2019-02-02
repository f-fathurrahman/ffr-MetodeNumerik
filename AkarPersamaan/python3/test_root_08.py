import numpy as np
from fixed_point import *

def f(x):
    return np.exp(-x)

x0 = 0.0
xroot, _ = fixed_point(f, x0, verbose=True)

print(xroot)
