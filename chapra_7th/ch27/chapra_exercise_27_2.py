import numpy as np
import matplotlib.pyplot as plt
from my_ode_solve import *
from my_root_solving import *

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

