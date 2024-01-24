from newton_raphson import *

def f(x):
    return x**3 - 35.0

def df(x):
    return 3*x**2

x, _ = newton_raphson(f, df, 3.0)

print("x = ", x)
print("x true = ", 35.0**(1/3))
print("x - xtrue = ", x - (35.0)**(1/3))