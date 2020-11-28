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
fxl = f(xl)
fxu = f(xu)
print("f(xl) = %f, f(xu) = %f" % (fxl, fxu))
xr = xu - fxu*(xl - xu)/(fxl - fxu)   # Eq. 5.7
fxr = f(xr)
print("xr = ", xr)
ε_t = abs(xr - x_true)/x_true*100 # error in percent
print("ε_t = %.1f%%" % ε_t)

if fxl*fxr < 0:
    xu = xr
else:
    xl = xr


xr_old = xr
# 2nd iteration
print("\n2nd iteration: ")
print("xl = %f, xu = %f" % (xl, xu))
fxl = f(xl)
fxu = f(xu)
print("f(xl) = %f, f(xu) = %f" % (fxl, fxu))
xr = xu - fxu*(xl - xu)/(fxl - fxu)   # Eq. 5.7
fxr = f(xr)
print("xr = ", xr)
ε_t = abs(xr - x_true)/x_true*100 # error in percent
print("ε_t = %.1f%%" % ε_t)
ε_a = abs(xr - xr_old)/xr*100 # error in percent
print("ε_a = %.1f%%" % ε_a)

if fxl*fxr < 0:
    xu = xr
else:
    xl = xr

xr_old = xr
# 3rd iteration
print("\n3rd iteration: ")
print("xl = %f, xu = %f" % (xl, xu))
fxl = f(xl)
fxu = f(xu)
print("f(xl) = %f, f(xu) = %f" % (fxl, fxu))
xr = xu - fxu*(xl - xu)/(fxl - fxu)   # Eq. 5.7
fxr = f(xr)
print("xr = ", xr)
ε_t = abs(xr - x_true)/x_true*100 # error in percent
print("ε_t = %.1f%%" % ε_t)
ε_a = abs(xr - xr_old)/xr*100 # error in percent
print("ε_a = %.1f%%" % ε_a)
