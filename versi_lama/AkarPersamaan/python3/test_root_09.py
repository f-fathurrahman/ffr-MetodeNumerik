import numpy as np
from newton_raphson import *
from ridder import *
from bisection import *

def f(x):
    return x**10 - 1

def df(x):
    return 10*x**9

x0 = 0.5

xroot, _ = newton_raphson(f, df, x0, verbose=True)

xroot, _ = bisection(f, 0.5, 1.25, verbose=True)

xroot, _ = ridder(f, 0.5, 1.25, verbose=True)

