from search_root_interval import *

def f(x):
    return x**3 - 10*x**2 + 5.0

x1 = 0.0
x2 = 1.0

for i in range(4):
    dx = (x2 - x1)/10.0
    x1, x2 = search_root_interval(f, x1, x2, dx)

x = 0.5*(x1 + x2)

print("x = {:6.4f}".format(x))