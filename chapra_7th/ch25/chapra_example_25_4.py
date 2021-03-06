import numpy as np
import matplotlib.pyplot as plt

def model_linear(t, v):
    g = 9.81
    c = 12.5
    m = 68.1
    return g - c*v/m

def model_nonlinear(t, v):
    g = 9.81
    c = 12.5
    m = 68.1
    a = 8.3
    b = 2.2
    vmax = 46
    return g - c/m*( v + a*(v/vmax)**b )

# One-step application of Euler's method for ODE
def ode_euler_1step(dfunc, xi, yi, h):
    return yi + dfunc(xi,yi)*h

def ode_euler(dfunc, x0, y0, h, Nstep):
    x = np.zeros(Nstep+1)
    y = np.zeros(Nstep+1)
    # Start with initial cond
    x[0] = x0
    y[0] = y0
    for i in range(0,Nstep):
        x[i+1] = x[i] + h
        y[i+1] = ode_euler_1step(dfunc, x[i], y[i], h)
    return x, y

def exact_sol(x):
    return -0.5*x**4 + 4*x**3 - 10*x**2 + 8.5*x + 1

# initial cond
t0 = 0.0
v0 = 0.0
tf = 15.0
h = 0.1
Nstep = int(tf/h)
t_linear, v_linear = ode_euler(model_linear, t0, v0, h, Nstep)
#
t_nonlinear, v_nonlinear = ode_euler(model_nonlinear, t0, v0, h, Nstep)


# Plot
plt.clf()
plt.plot(t_linear, v_linear, label="model linear")
plt.plot(t_nonlinear, v_nonlinear, label="model nonlinear")
plt.legend()
plt.savefig("IMG_example_25_4.png", dpi=150)
