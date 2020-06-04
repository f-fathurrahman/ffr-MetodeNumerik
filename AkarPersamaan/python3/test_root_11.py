from newton_raphson import *

def f(x):
    return x**4 - 38.5*x**3 + 458*x**2 - 1701.5*x + 741

def df(x):
    return 4*x**3 - 115.5*x**2 + 916*x - 1701.5

x, _ = newton_raphson(f, df, 12.0)
print("x = ", x)