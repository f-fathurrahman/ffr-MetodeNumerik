import numpy as np

m = 68.1 # mass, kg
v = 40.0 # velocity, m/s
t = 10.0 # time, s
g = 9.81

def f(c):
    return g*m/c * (1 - np.exp(-c/m*t)) - 40.0

x_true = 14.8011

# Initial guess
xl = 12.0
xu = 16.0

# First iteration
print("\n1st iteration: ")
print("xl = %f, xu = %f" % (xl, xu))
print("f(xl) = %f, f(xu) = %f" % (f(xl), f(xu)))
xr = (xl + xu)/2
print("xr = ", xr)
ε_t = abs(xr - x_true)/x_true*100 # error in percent
print("ε_t = %.1f %%" % ε_t)


if f(xl)*f(xr) < 0:
    xu = xr
else:
    xl = xr

# Second iteration
print("\n2nd iteration: ")
print("xl = %f, xu = %f" % (xl, xu))
print("f(xl) = %f, f(xu) = %f" % (f(xl), f(xu)))
xr = (xl + xu)/2
print("xr = ", xr)
ε_t = abs(xr - x_true)/x_true*100 # error in percent
print("ε_t = %.1f %%" % ε_t)


if f(xl)*f(xr) < 0:
    xu = xr
else:
    xl = xr

# Third iteration
print("\n3rd iteration: ")
print("xl = %f, xu = %f" % (xl, xu))
print("f(xl) = %f, f(xu) = %f" % (f(xl), f(xu)))
xr = (xl + xu)/2
print("xr = ", xr)
ε_t = abs(xr - x_true)/x_true*100 # error in percent
print("ε_t = %.1f %%" % ε_t)

