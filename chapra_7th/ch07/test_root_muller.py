from root_muller import *

def f(x):
    return x**3 - 13*x - 12.0

x1 = 1.1
xroot = root_muller(f, x1)
print("xroot = ", xroot)
if isinstance(xroot, complex):
    # Check for imaginary part
    if abs(xroot.imag) < 1e-12:
        # xroot is real
        xroot = xroot.real
