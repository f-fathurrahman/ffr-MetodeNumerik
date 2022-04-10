import numpy as np

def f(x):
    return np.sin(x)

a = 0.0
b = 1.0

Npoints = 100
xgrid = np.linspace(a, b, Npoints)
Δx = xgrid[1] - xgrid[0]
fgrid = f(xgrid)

I_exact = 1.0 - np.cos(1.0)
I = sum(fgrid)*Δx

print("Npoints = ", Npoints)
print("Integral approx = %18.10f" % I)
print("Error = %18.10e" % abs(I - I_exact))
