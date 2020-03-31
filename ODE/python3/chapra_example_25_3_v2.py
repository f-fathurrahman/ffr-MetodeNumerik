import numpy as np
import matplotlib.pyplot as plt

# the left hand side of dy/dx=... (in general depends on x and y)
# In the present case it only depends on x
def deriv(x, y):
    return -2*x**3 + 12*x**2 - 20*x + 8.5

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
x0 = 0.0
y0 = 1.0
xf = 4.0

# Using h=0.5
h = 0.5
Nstep = int(xf/h)
xs1, ys1 = ode_euler(deriv, x0, y0, h, Nstep)

# Using h=0.25
h = 0.25
Nstep = int(xf/h)
xs2, ys2 = ode_euler(deriv, x0, y0, h, Nstep)

# True solution
x_true = np.linspace(0.0, 4.0, 200)
y_true = exact_sol(x_true)

# Plot
plt.clf()
plt.plot(xs1, ys1, label="h=0.5", marker="o")
plt.plot(xs2, ys2, label="h=0.5", marker="o")
plt.plot(x_true, y_true, label="true") # Do not show the markers
plt.legend()
plt.savefig("IMG_example_25_3.png", dpi=150)
