import numpy as np
from math import ceil, log10
import matplotlib.pyplot as plt

# T == y[0]
# dT/dx == y[1]
# d2T/dx2 == dydx[1]
def deriv(x, y):
    #
    Nvec = len(y)
    assert Nvec == 2
    dydx = np.zeros(Nvec)
    #
    dydx[0] = y[1]
    dydx[1] = 0.15*y[0]
    #
    return dydx

def ode_rk4_1step(dfunc, xi, yi, h):
    k1 = dfunc(xi, yi)
    k2 = dfunc(xi + 0.5*h, yi + 0.5*k1*h)
    k3 = dfunc(xi + 0.5*h, yi + 0.5*k2*h)
    k4 = dfunc(xi + h, yi + k3*h)
    yip1 = yi + (k1 + 2*k2 + 2*k3 + k4)*h/6
    return yip1

def ode_solve(dfunc, do_1step, x0, y0, h, Nstep):
    Nvec = len(y0)
    x = np.zeros(Nstep+1)
    y = np.zeros((Nstep+1,Nvec))
    # Start with initial cond
    x[0] = x0
    y[0,:] = y0[:]
    for i in range(0,Nstep):
        x[i+1] = x[i] + h
        y[i+1,:] = do_1step(dfunc, x[i], y[i,:], h)
    return x, y

def obj_func(z0_guess):
    # Initial cond
    y0 = np.zeros(2)
    x0 = 0.0
    y0[0] = 240.0
    xf = 10.0  # end interval
    y0[1] = z0_guess
    Tf = 150.0
    #
    h = 0.5 # Step size
    Nstep = int( (xf-x0)/h )
    x, y = ode_solve(deriv, ode_rk4_1step, x0, y0, h, Nstep)
    # At the end of the interval
    Tf_guess = y[-1,0]
    return Tf_guess - Tf


def bisection(f, x1, x2, TOL=1.0e-9, NiterMax=None ):

    f1 = f(x1)
    if abs(f1) <= TOL:
        return x1, 0.0

    f2 = f(x2)
    if abs(f2) <= TOL:
        return x2, 0.0

    if f1*f2 > 0.0:
        raise RuntimeError("Root is not bracketed")

    # No NiterMax is provided
    # We calculate the default value here.
    if NiterMax == None:
        NiterMax = int(ceil( log10(abs(x2-x1)/TOL) )/ log10(2.0) ) + 10
        # extra 10 iterations

    # For the purpose of calculating relative error
    x3 = 0.0
    x3_old = 0.0

    print(13*" "+"Iter      Estimated          f(x)")
    print(13*" "+"----      ---------          ----")
    print("")

    for i in range(1,NiterMax+1):

        x3_old = x3
        x3 = 0.5*(x1 + x2)
        f3 = f(x3)

        print("bisection: %5d %18.10f %15.5e" % (i, x3, abs(f3)))

        if abs(f3) <= TOL:
            print("")
            print("bisection is converged in %d iterations" % i)
            # return the result
            return x3, abs(x3 - x3_old)

        if f2*f3 < 0.0:
            # sign of f2 and f3 is different
            # root is in [x2,x3]
            # change the interval bound of x1 to x3
            x1 = x3
            f1 = f3
        else:
            # sign of f1 and f3 is different
            # root is in [x1,x3]
            # change the interval bound of x2 to x3
            x2 = x3
            f2 = f3

    print("No root is found")
    return None, None


# Test plot the objective function
NpstPlt = 200
xplt = np.linspace(-500.0, 200.0, NpstPlt)
plt.clf()
yplt = np.zeros(NpstPlt)
for i in range(NpstPlt):
    yplt[i] = obj_func(xplt[i])
plt.plot(xplt, yplt)
plt.grid()
plt.savefig("IMG_exercise_27_2_obj_func.pdf")
#exit()

# For testing values of z0_1 and z0_2 which brackets obj_func
z0_1 = -200
Tf_1 = obj_func(z0_1)
z0_2 = 0.0
Tf_2 = obj_func(z0_2)
print("Tf_1 = ", Tf_1)
print("Tf_2 = ", Tf_2)

z0, Δ = bisection(obj_func, z0_1, z0_2, TOL=1.0e-9)

# Now solve the ODE with the obtained z0
y0 = np.zeros(2)
x0 = 0.0
y0[0] = 240.0
xf = 10.0  # end interval
y0[1] = z0
#
h = 0.5 # XXXX should be the same as the one used in obj_func
Nstep = int( (xf-x0)/h )
x, y = ode_solve(deriv, ode_rk4_1step, x0, y0, h, Nstep)
print("Tf = ", y[-1,0])

# Analytic solution, from exercise 27.1
def analytic_sols(x):
    return 3.01694393308884*np.exp(0.387298334620742*x) + \
    236.983056066911*np.exp(-0.387298334620742*x)


# Now plot the solution
plt.clf()
plt.plot(x, y[:,0], marker="o", label="Numeric")
plt.plot(x, analytic_sols(x), marker="x", label="Analytic")
plt.xlabel("x")
plt.ylabel("T")
plt.legend()
plt.savefig("IMG_exercise_27_2.pdf")

