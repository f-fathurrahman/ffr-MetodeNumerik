import numpy as np
import matplotlib.pyplot as plt
from math import exp

def ode_euler_1step(dfunc, xi, yi, h):
    return yi + dfunc(xi,yi)*h

def ode_heun_1step(dfunc, xi, yi, h):
    y0ip1 = yi + dfunc(xi,yi)*h
    avg = 0.5*( dfunc(xi,yi) + dfunc(xi+h,y0ip1) )*h
    return yi + avg

def ode_rk3_1step(dfunc, xi, yi, h):
    k1 = dfunc(xi, yi)
    k2 = dfunc(xi + 0.5*h, yi + 0.5*k1*h)
    k3 = dfunc(xi + h, yi - k1*h + 2*k2*h)
    yip1 = yi + (k1 + 4*k2 + k3)*h/6
    return yip1

def ode_rk4_1step(dfunc, xi, yi, h):
    k1 = dfunc(xi, yi)
    k2 = dfunc(xi + 0.5*h, yi + 0.5*k1*h)
    k3 = dfunc(xi + 0.5*h, yi + 0.5*k2*h)
    k4 = dfunc(xi + h, yi + k3*h)
    yip1 = yi + (k1 + 2*k2 + 2*k3 + k4)*h/6
    return yip1

def ode_rk5_1step(dfunc, xi, yi, h):
    k1 = dfunc(xi, yi)
    k2 = dfunc(xi + h/4, yi + k1*h/4)
    k3 = dfunc(xi + h/4, yi + k1*h/8 + k2*h/8)
    k4 = dfunc(xi + h/2, yi - k2*h/2 + k3*h)
    k5 = dfunc(xi + 3*h/4, yi + 3*k1*h/16 + 9*k4*h/16)
    k6 = dfunc(xi + h, yi - 3*k1*h/7 + 2*k2*h/7 + 12*k3*h/7 - 12*k4*h/7 + 8*k5*h/7)    
    yip1 = yi + (7*k1 + 32*k3 + 12*k4 + 32*k5 + 7*k6)*h/90
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

def deriv(x, y):
    return 4*exp(0.8*x) - 0.5*y

def exact_sol(x):
    return 4.0/1.3*( exp(0.8*x) - exp(-0.5*x) ) + 2*exp(-0.5*x)

# initial cond
x0 = 0.0
y0 = 2.0
# 
xf = 4.0

# These three arrays should have the same length
one_step_funcs = [ode_euler_1step, ode_heun_1step,
                  ode_rk3_1step, ode_rk4_1step, ode_rk5_1step]
labels = ["Euler", "Heun", "RK-3rd", "RK-4th", "Butcher"]
nf = [1, 2, 3, 4, 6]

Nmethods = len(one_step_funcs)

y_true = exact_sol(4.0)

# Loop is done over all methods, for a given g
# for debugging purpose
def do_solve_until_xf(h):
    Nstep = int(xf/h)
    ε_t = np.zeros(Nmethods)
    efforts = np.zeros(Nmethods)
    for i,f in enumerate(one_step_funcs):
        efforts[i] = nf[i]*(xf - x0)/h  # Equation (25.8.1)
        x, y = ode_solve(deriv, f, x0, y0, h, Nstep)
        ε_t[i] = abs(y_true - y[-1])/y_true * 100
        # Debug, make sure that x=4
        print("Method: %8s x = %f y = %f ε_t = %e%% effort = %8.1f" %
            (labels[i], x[-1], y[-1], ε_t[i], efforts[i]))


#do_solve_until_xf(2.0)
#do_solve_until_xf(1.0)
#do_solve_until_xf(0.5)
#do_solve_until_xf(0.25)
#do_solve_until_xf(0.125)
#do_solve_until_xf(0.1)


# Loop is done over all h (given as a list),
# for a given method and nf
def do_solve_until_xf(h, method, nf):
    Ndata = len(h)
    ε_t = np.zeros(Ndata)
    efforts = np.zeros(Ndata)
    for i in range(Ndata):
        Nstep = int(xf/h[i])
        efforts[i] = nf*(xf - x0)/h[i]  # Equation (25.8.1)
        x, y = ode_solve(deriv, method, x0, y0, h[i], Nstep)
        ε_t[i] = abs(y_true - y[-1])/y_true * 100
    return efforts, ε_t

h = [2.0, 1.0, 0.5, 0.25, 0.125, 0.1]
plt.clf()
for i,f in enumerate(one_step_funcs):
    efforts, ε_t = do_solve_until_xf(h, f, nf[i])
    plt.plot(efforts, np.log10(ε_t), marker="o", label=labels[i])
plt.xlabel("Effort")
plt.ylabel("log10(ε_t)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("IMG_chapra_example_25_8.pdf")
