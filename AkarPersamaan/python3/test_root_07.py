from newton_raphson import *

def f(x):
    return x**4 - 6.4*x**3 + 6.45*x**2 + 20.538*x - 31.752

def df(x):
    return 4.0*x**3 - 19.2*x**2 + 12.9*x + 20.538

x0 = 3.0 # initial guess of root
x1, _ = newton_raphson(f, df, x0, verbose=True, TOL=1e-12, multiplicity=2)

print("x1 = ", x1)
print("fx1 = ", f(x1))

x0 = 2.0 # initial guess of root
x2, _ = newton_raphson(f, df, x0, verbose=True, TOL=1e-12)

print("x2 = ", x2)
print("fx2 = ", f(x2))

print("diff = ", x2-x1)