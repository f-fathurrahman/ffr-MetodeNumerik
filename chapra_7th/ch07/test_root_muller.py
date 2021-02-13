from cmath import sin, cos
from root_muller import *

def f(x):
    return x**3 - 13*x - 12.0

#x1 = 1.0
#xroot = root_muller(f, x1)

# Chapra exercise 6.7
def chapra_exe_6_7(x):
    return sin(x) + cos(1 + x**2) - 1

x1 = 1.0 # complex root?
#x1 = 2.0
#x1 = 10.0
xr = root_muller(chapra_exe_6_7, x1)