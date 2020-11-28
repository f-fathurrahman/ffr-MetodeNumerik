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

NiterMax = 10

xr_old = 0.0
ε_a = 100.0

for iiter in range(1,NiterMax+1):
    #
    xr = (xl + xu)/2
    #
    print("%3d" % iiter, end="")
    print("%10.4f %10.4f" % (xl, xu), end="")
    print("%10.4f" % xr, end="")
    #
    if iiter > 1:
        ε_a = abs( (xr - xr_old) / xr )*100 # in percent
        print("%10.3f%%" % ε_a, end="")
    else: # to make the print out looks tidy
        print("           ", end="")
    #
    ε_t = abs( (xr - x_true) / x_true )*100 # in percent
    print("%10.3f%%" % ε_t)
    #
    if ε_a < 0.5:
        print("Convergence achieved: ε_a < 0.5%")
        break
    #
    if f(xl)*f(xr) < 0:
        xu = xr
    else:
        xl = xr
    #
    xr_old = xr

