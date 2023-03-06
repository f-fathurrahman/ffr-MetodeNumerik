import numpy as np
from scipy.optimize import minimize

def my_func(x):
    return 2*np.sin(x) - x**2/10
    # The equation shows x**2/2
    # However the MATLAB code uses x**2/10
    # The result is shown using x**2/10

def obj_func(x):
    return -my_func(x)

x0 = 2.0 # initial search point, somewhat arbitrary
# here we use the midpoint between x=0 and xu=4 that are given in the book

result = minimize(obj_func, x0, method="Nelder-Mead") # default method is BFGS

# Full output (might differ for each methods)
print("result = ")
print(result)

# The result
print("result.x (arg) = ", result.x)
print("result.fun (max value) = ", -result.fun)
