import numpy as np

def deriv(x, y):

    a = 100.0
    b = 15.0

    Nvec = len(y)
    assert Nvec == 2

    # Output array
    dydx = np.zeros(Nvec)
    
    # remember that in Python the array starts at 0
    # y1 = y[0] = θ
    # y2 = y[1] = θdot
    dydx[0] = y[1]
    dydx[1] = ( a*(b - y[0]) - y[0]*y[1]**2 ) / (1 + y[0]**2)
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
t0 = 0.0
theta0 = np.array([2*np.pi, 0])

t = 0.5
h = 0.01

Nstep = int(t/h)

t, theta = ode_solve(deriv, ode_rk4_1step, t0, theta0, h, Nstep)
print("")
print("--------------------------------")
print(" t         theta1         theta2")
print("--------------------------------")
for i in range(len(t)):
    print("%5.1f %10.6f %10.6f" % (t[i], theta[i,0], theta[i,1]))

import matplotlib.pyplot as plt
plt.plot(t, theta[:,0], label="theta", marker="o")
plt.plot(t, theta[:,1], label="thetadot", marker="o")
plt.legend()
plt.grid()
plt.savefig("IMG_Tes_2020_Soal_5.pdf")
