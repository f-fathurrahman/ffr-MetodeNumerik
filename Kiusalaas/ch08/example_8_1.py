import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.size": 14
})

from my_ode_solve import *
from my_root_solving import *

def deriv(x, y):
    #
    Nvec = len(y)
    assert Nvec == 2
    dydx = np.zeros(Nvec)
    dydx[0] = y[1]
    dydx[1] = -3*y[0]*y[1]
    #
    return dydx


def obj_func(z0_guess):
    # Initial cond
    x0 = 0.0
    y0 = np.zeros(2)
    #
    y0[0] = 0.0
    y0[1] = z0_guess
    #
    xf = 2.0
    yf = 1.0
    #
    h = 0.1 # Step size
    Nstep = int( (xf-x0)/h )
    x, y = ode_solve(deriv, ode_rk4_1step, x0, y0, h, Nstep)
    # At the end of the interval
    yf_guess = y[-1,0]
    return yf_guess - yf


# For testing values of z0_1 and z0_2 which brackets obj_func
z0_1 = 1.0
y_1 = obj_func(z0_1)
z0_2 = 2.0
y_2 = obj_func(z0_2)
print("y_1 = ", y_1)
print("y_2 = ", y_2)


z0 = root_bisection(obj_func, z0_1, z0_2, TOL=1.0e-9)
#z0 = root_regula_falsi(obj_func, z0_1, z0_2, TOL=1.0e-9)

# Now solve the ODE with the obtained z0
x0 = 0.0
y0 = np.zeros(2)
y0[0] = 0.0
y0[1] = z0
xf = 2.0
#
h = 0.1 # Step size
Nstep = int( (xf-x0)/h )
x, y = ode_solve(deriv, ode_rk4_1step, x0, y0, h, Nstep)
print("yf = ", y[-1,0])  # CHECK: Should give a value close to 200.0

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

