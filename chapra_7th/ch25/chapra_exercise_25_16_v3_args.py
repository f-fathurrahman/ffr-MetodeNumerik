from math import sin
import numpy as np
import matplotlib.pyplot as plt

def ode_rk4_1step(dfunc, xi, yi, h, args=()):
    k1 = dfunc(xi, yi, *args)
    k2 = dfunc(xi + 0.5*h, yi + 0.5*k1*h, *args)
    k3 = dfunc(xi + 0.5*h, yi + 0.5*k2*h, *args)
    k4 = dfunc(xi + h, yi + k3*h, *args)
    yip1 = yi + (k1 + 2*k2 + 2*k3 + k4)*h/6
    return yip1

def ode_solve(dfunc, do_1step, x0, y0, h, Nstep, args=()):
    Nvec = len(y0)
    x = np.zeros(Nstep+1)
    y = np.zeros((Nstep+1,Nvec))
    # Start with initial cond
    x[0] = x0
    y[0,:] = y0[:]
    for i in range(0,Nstep):
        x[i+1] = x[i] + h
        y[i+1,:] = do_1step(dfunc, x[i], y[i,:], h, args=args)
    return x, y


def deriv(x, y, m, k, c):
    Nvec = len(y)
    # Here we make an assertion to make sure that y is a 2-component vector
    # Uncomment this line if the code appears to be slow
    assert Nvec == 2
    # Output array
    dydx = np.zeros(Nvec)
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


m = 0.1
k = 40.0
c_under = 0.5
c_cri   = 4.0
c_over  = 8.0

# ordering of m,k,c is important (refer to the definition of deriv)

t, x_under = ode_solve(deriv, ode_rk4_1step, t0, x0, Δt, Nstep, args=(m,k,c_under))
t, x_cri   = ode_solve(deriv, ode_rk4_1step, t0, x0, Δt, Nstep, args=(m,k,c_cri))
t, x_over  = ode_solve(deriv, ode_rk4_1step, t0, x0, Δt, Nstep, args=(m,k,c_over))

plt.clf()
plt.plot(t, x_under[:,0], label="underdamped")
plt.plot(t, x_cri[:,0], label="critical")
plt.plot(t, x_over[:,0], label="overdamped")
plt.grid(True)
plt.xlabel("t (s)")
plt.ylabel("Displacement (m)")
plt.legend()
plt.tight_layout()
plt.savefig("IMG_chapra_exercise_25_16_v3_args.pdf")
