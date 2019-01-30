import math
from search_root_interval import *
from bisection import *

def f(x):
    return x - math.tan(x)

a = 0.0
b = 20.0
dx = 0.01

while True:
    x1, x2 = search_root_interval(f,a,b,dx)
    if x1 != None:
        a = x2
        root, _ = bisection(f, x1, x2)
        if root != None:
            print("Root is found: {:18.10f}".format(root))
    else:
        print("\nDone")
        break

