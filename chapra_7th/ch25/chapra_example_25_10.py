import numpy as np

def deriv(x, y):
    Nvec = len(y)
    # Here we make an assertion to make sure that y is a 2-component vector
    # Uncomment this line if the code appears to be slow
    assert Nvec == 2
    # Output array
    dydx = np.zeros(Nvec)
    # remember that in Python the array starts at 0
    # y1 = y[0]
    # y2 = y[1]
    dydx[0] = -0.5*y[0]
    dydx[1] = 4 - 0.3*y[1] - 0.1*y[0]
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
y0 = np.array([4.0, 6.0])
h = 0.5
Nstep = 4
x, y = ode_solve(deriv, ode_rk4_1step, x0, y0, h, Nstep)
print("")
print("---------------------------")
print(" x         y1         y2")
print("---------------------------")
for i in range(len(x)):
    print("%5.1f %10.6f %10.6f" % (x[i], y[i,0], y[i,1]))
