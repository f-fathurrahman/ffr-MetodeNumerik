import numpy as np

def f(x):
    return np.sin(x)

a = 0.0
b = 1.0

I_exact = 1.0 - np.cos(1.0)
for Npoints in [10, 100, 1000, 5000, 10000]:
    xgrid = np.linspace(a, b, Npoints)
    Δx = xgrid[1] - xgrid[0]
    fgrid = f(xgrid)
    I = sum(fgrid)*Δx
    error = abs(I - I_exact)
    print("%8d %15.10f %10.5e" % (Npoints, I, error))
