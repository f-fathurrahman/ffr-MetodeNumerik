import numpy as np
from golden_section_optim import *

def f(x):
    return 2*np.sin(x) - x**2/10

xopt, fx = golden_section_optim(f, 0.0, 4.0)
print("xopt    = ", xopt)
print("f(xopt) = ", fx)