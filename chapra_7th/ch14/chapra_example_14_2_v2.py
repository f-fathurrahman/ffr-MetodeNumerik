import numpy as np
import matplotlib.pyplot as plt

def my_func(x, y):
    return 2*x*y + 2*x - x**2 - 2*y**2

def grad_my_func(x, y):
    dfdx = 2*y + 2 - 2*x
    dfdy = 2*x - 4*y
    return np.array([dfdx, dfdy]) # return as numpy array

# Calc the gradient
x0, y0 = -1.0, 1.0
gx, gy = grad_my_func(x0, y0)


Δ = 0.1

# Now suppose that we want to search for MAXIMUM value
# We want to search for new point: (x0 + Δ*gx, y0 + Δ*gy)
# or we move in the direction of STEEPEST ASCENT
# which gives larger value than current value

# Value of the function at (x0,y0)
f1 = my_func(x0, y0)
f2 = my_func(x0 + Δ*gx, y0 + Δ*gy)
print("\nTrying to find maximum")
print("Old value: %18.10f" % f1)
print("New value: %18.10f" % f2)
if f2 > f1:
    print("Good: Function value is increasing")
else:
    print("Bad: Function value is decreasing")
    print("Step length is too large") # or Δ is too large




# Now suppose that we want to search for MINIMUM value
# We want to search for new point: (x0 - Δ*gx, y0 - Δ*gy)
# or we move in the direction of STEPEEST DESCENT
# which gives smaller value than current value

# Value of the function at (x0,y0)
f1 = my_func(x0, y0)
f2 = my_func(x0 - Δ*gx, y0 - Δ*gy)
print("\nTrying to find minimum")
print("Old value: %18.10f" % f1)
print("New value: %18.10f" % f2)
if f2 < f1:
    print("Good: Function value is decreasing")
else:
    print("Bad: Function value is increasing")
    print("Step length is too large") # or Δ is too large
