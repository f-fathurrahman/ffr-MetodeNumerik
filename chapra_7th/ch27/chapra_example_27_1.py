# Shooting method for linear BVP

# d2T/dx2 + h'(T_a - T) = 0
# Boundary condition:
#   T(0) = 40
#   T(10) = 200

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# T == y[0]
# dT/dx == y[1]
# d2T/dx2 == dydx[1]
def deriv(x, y):
    #
    Nvec = len(y)
    assert Nvec == 2
    dydx = np.zeros(Nvec)
    #
    h = 0.01
    T_a = 20.0
    #
    dydx[0] = y[1]
    dydx[1] = h*(y[0] - T_a)
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

# Initial cond
x0 = 0.0
y0 = np.zeros(2)
y0[0] = 40.0 # from the boundary condition, T(0) = 40
z0_1 = 10.0 # Guess for z = dT/dx == y[1], save to it variable for later use
y0[1] = z0_1
xf = 10.0  # end interval
Tf = 200.0 # Boundary condition, T(10) = 200

h = 2.0 # Step size
Nstep = int( (xf-x0)/h )
x, y = ode_solve(deriv, ode_rk4_1step, x0, y0, h, Nstep)
# At the end of the interval
Tf_1 = y[-1,0]
print("First guess: T(10) = y[-1,0] = ", Tf_1)

# Integrate again, now with new guess for z(0) = y0[1]
z0_2 = 20.0
y0[1] = z0_2 # Guess for z = dT/dx == y[1]
x, y = ode_solve(deriv, ode_rk4_1step, x0, y0, h, Nstep)
# At the end of the interval
Tf_2 = y[-1,0]
print("Second guess: T(10) = y[-1,0] = ", Tf_2)

# Using linear interp to guess what value of z(0) which gives T(10) = 200
z0_new = z0_1 + (z0_2 - z0_1)/(Tf_2 - Tf_1) * (Tf - Tf_1)
print("z0_new = ", z0_new)

# Now solve the IVP using z0_new
y0[1] = z0_new # Guess for z = dT/dx == y[1]
x, y = ode_solve(deriv, ode_rk4_1step, x0, y0, h, Nstep)
# At the end of the interval
Tf_3 = y[-1,0]
print("Third guess: T(10) = y[-1,0] = ", Tf_3) # Should give Tf = 200

# The obtained solution, i.e. y, is the solution to the BVP
# Let's visualize the solution
plt.clf()
plt.plot(x, y[:,0], marker="o", label="Temperature")
plt.xlabel("x")
plt.ylabel("T")
plt.legend()
plt.savefig("IMG_example_27_1.png", dpi=150)

# Of course, you may try to visualize the y's that's produced
# by guessing z0 for comparison.
