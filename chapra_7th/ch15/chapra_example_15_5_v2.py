import numpy as np
from scipy.optimize import minimize_scalar

def my_func(x):
    return 2*np.sin(x) - x**2/10
    # The equation shows x**2/2
    # However the MATLAB code uses x**2/10
    # The result is shown using x**2/10

def obj_func(x):
    return -my_func(x)

xl = 0.0
xu = 4.0

result = minimize_scalar(obj_func, bounds=(xl,xu))

# Full output (might differ for each methods)
print("result = ")
print(result)

# The result
print("result.x (arg) = ", result.x)
print("result.fun (max value) = ", -result.fun)
