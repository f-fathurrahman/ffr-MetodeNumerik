import numpy as np

def deriv(x, y):
    
    m = 80.0
    Cd = 0.2028
    g = 9.8

    Nvec = len(y)
    assert Nvec == 2
    dydx = np.zeros(Nvec)

    dydx[0] = y[1]
    dydx[1] = -g + Cd/m * y[1]**2
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
y0 = np.array([500.0, 0])

t = 13.0
h = 0.1

Nstep = int(t/h)

t, y = ode_solve(deriv, ode_rk4_1step, t0, y0, h, Nstep)
print("")
print("--------------------------------")
print(" t         theta1         theta2")
print("--------------------------------")
for i in range(len(t)):
    print("%5.1f %10.6f %10.6f" % (t[i], y[i,0], y[i,1]))

import matplotlib.pyplot as plt
plt.plot(t, y[:,0], label="y", marker="o")
plt.plot(t, y[:,1], label="vel", marker="o")
plt.legend()
plt.grid()
plt.savefig("IMG_Tes_2020_Soal_5_jatuh.pdf")
