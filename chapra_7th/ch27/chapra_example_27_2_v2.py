# Shooting method for nonlinear BVP

# d2T/dx2 + h'(T_a - T)**4 = 0
# Boundary condition:
#   T(0) = 40
#   T(10) = 200

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.size": 14
})

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
    h = 5e-8
    T_a = 20.0
    #
    dydx[0] = y[1]
    dydx[1] = h*(y[0] - T_a)**4
    #
    return dydx


def obj_func(z0_guess):
    # Initial cond
    x0 = 0.0
    y0 = np.zeros(2)
    y0[0] = 40.0 # from the boundary condition, T(0) = 40
    y0[1] = z0_guess
    xf = 10.0  # end interval
    Tf = 200.0 # Boundary condition, T(10) = 200
    #
    h = 2.0 # Step size
    Nstep = int( (xf-x0)/h )
    x, y = ode_solve(deriv, ode_rk4_1step, x0, y0, h, Nstep)
    # At the end of the interval
    Tf_guess = y[-1,0]
    return Tf_guess - Tf


# For testing values of z0_1 and z0_2 which brackets obj_func
z0_1 = 5.0
Tf_1 = obj_func(z0_1)
z0_2 = 11.0
Tf_2 = obj_func(z0_2)
print("Tf_1 = ", Tf_1)
print("Tf_2 = ", Tf_2)

z0 = root_bisection(obj_func, z0_1, z0_2, TOL=1.0e-9)
#z0 = root_regula_falsi(obj_func, z0_1, z0_2, TOL=1.0e-9)

# Now solve the ODE with the obtained z0
x0 = 0.0
y0 = np.zeros(2)
y0[0] = 40.0 # from the boundary condition, T(0) = 40
y0[1] = z0
xf = 10.0  # end interval
#
h = 2.0 # Step size
Nstep = int( (xf-x0)/h )
x, y = ode_solve(deriv, ode_rk4_1step, x0, y0, h, Nstep)
print("Tf = ", y[-1,0])  # CHECK: Should give a value close to 200.0

# Now plot the solution
plt.clf()
plt.plot(x, y[:,0], marker="o", label="Temperature")
plt.xlabel("x")
plt.ylabel("T")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("IMG_example_27_2.pdf")
plt.savefig("IMG_example_27_2.png", dpi=150)

