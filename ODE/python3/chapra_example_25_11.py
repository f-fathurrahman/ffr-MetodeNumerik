from math import sin
import numpy as np
import matplotlib.pyplot as plt

def pendulum_ode(x, y):
    Nvec = len(y)
    # Here we make an assertion to make sure that y is a 2-component vector
    # Uncomment this line if the code appears to be slow
    assert Nvec == 4
    # Output array
    dydx = np.zeros(Nvec)
    # remember that in Python the array starts at 0
    # y1 = y[0]
    # y2 = y[1], etc ...
    dydx[0] = y[1]
    dydx[1] = -16.1*y[0]
    # Nonlinear effect
    dydx[2] = y[3]
    dydx[3] = -16.1*sin(y[2])
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

# initial cond
x0 = 0.0
# Small displacement
y0 = np.array([0.1, 0.0, 0.1, 0.0])
h = 0.01 # try playing with this parameter
xf = 4.0
Nstep = int(xf/h)
x, y = ode_solve(pendulum_ode, ode_rk4_1step, x0, y0, h, Nstep)

plt.clf()
plt.plot(x, y[:,0], label="y1")
plt.plot(x, y[:,1], label="y2")
plt.plot(x, y[:,2], label="y3")
plt.plot(x, y[:,3], label="y4")
plt.title("Small displacement case")
plt.ylim(-4,4) # The same for both small and large displacement (to match the figure in the book)
plt.legend()
plt.savefig("IMG_example_25_11_v1.png", dpi=150)


# initial cond
x0 = 0.0
# Large displacement
y0 = np.array([np.pi/4, 0.0, np.pi/4, 0.0])
h = 0.01 # try playing with this parameter
xf = 4.0
Nstep = int(xf/h)
x, y = ode_solve(pendulum_ode, ode_rk4_1step, x0, y0, h, Nstep)

plt.clf()
plt.plot(x, y[:,0], label="y1")
plt.plot(x, y[:,1], label="y2")
plt.plot(x, y[:,2], label="y3")
plt.plot(x, y[:,3], label="y4")
plt.title("Large displacement case")
plt.ylim(-4,4) # The same for both small and large displacement (to match the figure in the book)
plt.legend()
plt.savefig("IMG_example_25_11_v2.png", dpi=150)
