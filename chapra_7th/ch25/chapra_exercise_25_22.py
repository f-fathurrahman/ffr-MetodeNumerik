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

def dfunc_for_v(x, y):
    g0 = 9.81
    R = 6.37e6
    return -g0*R**2/(R + x)**2

t0 = 0.0
v0 = 1500.0

Δt = 0.05
Nstep = int(200.0/Δt)
t, v = ode_solve(dfunc_for_v, ode_rk4_1step, t0, [v0], Δt, Nstep)

plt.clf()
plt.plot(t, v)
plt.grid(True)
plt.savefig("IMG_chapra_exercise_22.pdf")