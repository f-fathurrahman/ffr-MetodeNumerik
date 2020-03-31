import numpy as np
import matplotlib.pyplot as plt

def deriv(x, y):
    return -2*x**3 + 12*x**2 - 20*x + 8.5

#
def ode_euler_1step(dfunc, xi, yi, h):
    return yi + dfunc(xi,yi)*h

def ode_heun_1step(dfunc, xi, yi, h):
    y0ip1 = yi + dfunc(xi,yi)*h
    avg = 0.5*( dfunc(xi,yi) + dfunc(xi+h,y0ip1) )*h
    return yi + avg

def ode_midpoint_1step(dfunc, xi, yi, h):
    yip12 = yi + deriv(xi,yi)*h/2  # midpoint value
    xip12 = xi + 0.5*h             # midpoint
    yip1 = yi + deriv(xip12,yip12)*h
    return yip1

def ode_ralston_1step(dfunc, xi, yi, h):
    k1 = dfunc(xi, yi)
    k2 = dfunc(xi + 3*h/4, yi + 3*k1*h/4)
    yip1 = yi + (k1/3 + 2*k2/3)*h
    return yip1

# do_1step is a function, in the present case it
# may be one of 1-step ode solver defined above
def ode_solve(dfunc, do_1step, x0, y0, h, Nstep):
    x = np.zeros(Nstep+1)
    y = np.zeros(Nstep+1)
    # Start with initial cond
    x[0] = x0
    y[0] = y0
    for i in range(0,Nstep):
        x[i+1] = x[i] + h
        y[i+1] = do_1step(dfunc, x[i], y[i], h)
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

# You can pass function as an argument
# In the following, we pass 1step ode solvers (functions) to ode_solve
# You can even use loop to simplify the codes below

# Euler
x_euler, y_euler = ode_solve(deriv, ode_euler_1step, x0, y0, h, Nstep)

# Heun
x_heun, y_heun = ode_solve(deriv, ode_heun_1step, x0, y0, h, Nstep)

# Midpoint
x_midpoint, y_midpoint = ode_solve(deriv, ode_midpoint_1step, x0, y0, h, Nstep)

# Ralston
x_ralston, y_ralston = ode_solve(deriv, ode_ralston_1step, x0, y0, h, Nstep)

# True solution
x_true = np.linspace(0.0, 4.0, 200)
y_true = exact_sol(x_true)

# Plot
plt.clf()
plt.plot(x_euler, y_euler, label="Euler", marker="o")
plt.plot(x_heun, y_heun, label="Heun", marker="o")
plt.plot(x_midpoint, y_midpoint, label="Midpoint", marker="o")
plt.plot(x_ralston, y_ralston, label="Ralston", marker="o")
plt.plot(x_true, y_true, label="true") # Do not show the markers
plt.legend()
plt.savefig("IMG_example_25_6.png", dpi=150)


# Table stuffs
x = x_heun # the x are the same for all methods, we use the Heun's x
Nx = len(x)
# print header
print("")
print("----------------------------------------------------------------------")
print("                      Heun              Midpoint           Ralston    ")
print(" x    y_true       y       |ε_t|      y       |ε_t|      y       |ε_t|")
print("----------------------------------------------------------------------")
for i in range(Nx):
    #
    y_true = exact_sol(x[i])
    #
    err_heun = abs(y_true - y_heun[i])/y_true * 100
    err_midpoint = abs(y_true - y_midpoint[i])/y_true * 100
    err_ralston = abs(y_true - y_ralston[i])/y_true * 100
    #
    print("%.1f   %.5f   %.5f   %5.2f%%   %.5f   %5.2f%%   %.5f   %5.2f%%" % (
        x[i], y_true, y_heun[i], err_heun, 
        y_midpoint[i], err_midpoint,
        y_ralston[i], err_ralston) )
