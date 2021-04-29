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
    C1 = 1e-7
    y0s = y[0] + 273
    # y[0] -> y0s - 273
    Tb = 150.0 # Linearization energy
    linear_term = C1*(Tb+273.0)**4 + 4*C1*(Tb+273.0)**3*(y[0] - Tb)
    dydx[1] = linear_term - 4*(150 - y[0])
    #
    return dydx

def obj_func(z0_guess):
    # Initial cond
    y0 = np.zeros(2)
    x0 = 0.0
    y0[0] = 200.0
    xf = 0.5  # end interval
    y0[1] = z0_guess
    Tf = 100.0
    #
    h = 0.05 # Step size
    Nstep = int( (xf-x0)/h )
    x, y = ode_solve(deriv, ode_rk4_1step, x0, y0, h, Nstep)
    # At the end of the interval
    Tf_guess = y[-1,0]
    #print("Tf_guess = ", Tf_guess)
    return Tf_guess - Tf


# Test plot the objective function
NpstPlt = 200
xplt = np.linspace(-1000.0, -400.0, NpstPlt)
plt.clf()
yplt = np.zeros(NpstPlt)
for i in range(NpstPlt):
    yplt[i] = obj_func(xplt[i])
plt.plot(xplt, yplt)
plt.grid()
plt.savefig("IMG_exercise_27_7_obj_func.pdf")
#exit()

# For testing values of z0_1 and z0_2 which brackets obj_func
z0_1 = -900.0
Tf_1 = obj_func(z0_1)
z0_2 = -800.0
Tf_2 = obj_func(z0_2)
print("Tf_1 = ", Tf_1)
print("Tf_2 = ", Tf_2)
#exit()

z0 = root_bisection(obj_func, z0_1, z0_2, TOL=1.0e-9, NiterMax=100)
#z0 = root_regula_falsi(obj_func, z0_1, z0_2, TOL=1.0e-9, NiterMax=100)

# Now solve the ODE with the obtained z0
y0 = np.zeros(2)
x0 = 0.0
y0[0] = 200.0
xf = 0.5  # end interval
y0[1] = z0
#
h = 0.05 # XXXX should be the same as the one used in obj_func
Nstep = int( (xf-x0)/h )
x, y = ode_solve(deriv, ode_rk4_1step, x0, y0, h, Nstep)
print("Tf = ", y[-1,0])

# Now plot the solution
plt.clf()
plt.plot(x, y[:,0], marker="o", label="Numeric")
#plt.plot(x, analytic_sols(x), marker="x", label="Analytic")
plt.xlabel("x")
plt.ylabel("y")
plt.ylim(80,210)
plt.legend()
plt.grid(True)
plt.title("Linear PDE")
plt.savefig("IMG_exercise_27_7.pdf")

