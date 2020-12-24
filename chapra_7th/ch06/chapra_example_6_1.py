# Simple fixed-point iteration

# f(x) = exp(-x) - x
# The equation is transformed to
# x = exp(-x) = g(x)
#
# x_{i+1} = exp(-x_{i})

import numpy as np

def g(x):
    return np.exp(-x)

# Initial guess
x = 0.0
x_true = 0.56714329

print()
print("Initial point: x = ", x)
print()
for i in range(1,11):
    xnew = g(x)
    ε_a = np.abs( (xnew - x)/xnew )*100 # in percent
    ε_t = np.abs( (x_true - xnew)/x_true )*100
    print("%3d %10.6f %10.2f%% %10.2f%%" % (i, xnew, ε_a, ε_t))
    x = xnew


