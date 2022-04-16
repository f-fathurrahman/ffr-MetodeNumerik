from math import sin
import numpy as np
import matplotlib.pyplot as plt

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


def deriv_underdamped(x, y):
    # Parameters (local)
    m = 0.1   # kg
    k = 40.0   # N/m
    c = 0.5    # Ns/m underdamped
    Nvec = len(y)
    # Here we make an assertion to make sure that y is a 2-component vector
    # Uncomment this line if the code appears to be slow
    assert Nvec == 2
    # Output array
    dydx = np.zeros(Nvec)
    # c, k, m are accessed as global parameters
    dydx[0] = y[1]
    dydx[1] = -(c*y[1] + k*y[0])/m
    #
    return dydx


def deriv_critical_damped(x, y):
    # Parameters (local)
    m = 0.1   # kg
    k = 40.0   # N/m
    c = 4.0    # Ns/m critically damped
    Nvec = len(y)
    # Here we make an assertion to make sure that y is a 2-component vector
    # Uncomment this line if the code appears to be slow
    assert Nvec == 2
    # Output array
    dydx = np.zeros(Nvec)
    # c, k, m are accessed as global parameters
    dydx[0] = y[1]
    dydx[1] = -(c*y[1] + k*y[0])/m
    #
    return dydx


def deriv_overdamped(x, y):
    # Parameters (local)
    m = 0.1   # kg
    k = 40.0   # N/m
    c = 8.0    # Ns/m overdamped
    Nvec = len(y)
    # Here we make an assertion to make sure that y is a 2-component vector
    # Uncomment this line if the code appears to be slow
    assert Nvec == 2
    # Output array
    dydx = np.zeros(Nvec)
    # c, k, m are accessed as global parameters
    dydx[0] = y[1]
    dydx[1] = -(c*y[1] + k*y[0])/m
    #
    return dydx

# Below, we will use:
# - t as independent variable
# - x as dependent variable or the function we are trying to find
#   by solving the ODE

# initial conditions
t0 = 0.0
initial_displ = 0.1
initial_vel = 10.0
x0 = np.array([initial_displ, initial_vel])
Δt = 0.001 # try playing with this parameter
tf = 2.0
Nstep = int(tf/Δt)
print("Nstep = ", Nstep)

t, x_underdamped = ode_solve(deriv_underdamped, ode_rk4_1step, t0, x0, Δt, Nstep)
t, x_critical_damped = ode_solve(deriv_critical_damped, ode_rk4_1step, t0, x0, Δt, Nstep)
t, x_overdamped = ode_solve(deriv_overdamped, ode_rk4_1step, t0, x0, Δt, Nstep)

plt.clf()
plt.plot(t, x_underdamped[:,0], label="underdamped")
plt.plot(t, x_critical_damped[:,0], label="critical")
plt.plot(t, x_overdamped[:,0], label="overdamped")
plt.grid(True)
plt.xlabel("t (s)")
plt.ylabel("Displacement (m)")
plt.legend()
plt.tight_layout()
plt.savefig("IMG_chapra_exercise_25_16_v3.pdf")
