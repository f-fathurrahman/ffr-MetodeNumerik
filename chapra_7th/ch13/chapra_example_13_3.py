import numpy as np

def my_func(x):
    return 2*np.sin(x) - x**2/10

def d1_my_func(x):
    return 2*np.cos(x) - x/5

def d2_my_func(x):
    return -2*np.sin(x) - 1/5

SMALL = np.finfo(np.float64).resolution
NiterMax = 100

# Initial guess
x0 = 2.5
fopt_old = np.nan

for iiter in range(1,NiterMax+1):
    f0 = my_func(x0)
    df0 = d1_my_func(x0)
    d2f0 = d2_my_func(x0)
    if abs(df0) > SMALL:
        x1 = x0 - df0/d2f0
        f1 = my_func(x1)
        print("%18.10f %18.10f %18.10e" % (x1, f1, abs(f1 - fopt_old)))
    else:
        print("Converged")
        break
    x0 = x1
    fopt_old = f1

