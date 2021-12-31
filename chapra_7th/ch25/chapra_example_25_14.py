import numpy as np

def ode_rk4_1step(dfunc, xi, yi, h):
    k1 = dfunc(xi, yi)
    k2 = dfunc(xi + 0.5*h, yi + 0.5*k1*h)
    k3 = dfunc(xi + 0.5*h, yi + 0.5*k2*h)
    k4 = dfunc(xi + h, yi + k3*h)
    yip1 = yi + (k1 + 2*k2 + 2*k3 + k4)*h/6
    return yip1

def ode_rk4(dfunc, x0, y0, h, Nstep):
    x = np.zeros(Nstep+1)
    y = np.zeros(Nstep+1)
    #
    x[0] = x0
    y[0] = y0
    for i in range(0,Nstep):
        x[i+1] = x[i] + h
        y[i+1] = ode_rk4_1step(dfunc, x[i], y[i], h)
    return x, y

def ode_rk45_1step(dfunc, xi, yi, h):
    k1 = dfunc(xi, yi)
    k2 = dfunc(xi + h/5, yi + k1*h/5)
    k3 = dfunc(xi + 3*h/10, yi + 3*k1*h/40 + 9*k2*h/40)
    k4 = dfunc(xi + 3*h/5, yi + 3*k1*h/10 - 9*k2*h/10 + 6*k3*h/5)
    k5 = dfunc(xi + h, yi - 11*k1*h/54 + 5*k2*h/2 - 70*k3*h/27 + 35*k4*h/27)
    k6 = dfunc(xi + 7*h/8, yi + 1631*k1*h/55296 + 175*k2*h/512 + 575*k3*h/13824 +
              44275*k4*h/110592 + 253*k5*h/4096)
    yip1 = yi + (37*k1/378 + 250*k3/621 + 125*k4/594 + 512*k6/1771)*h
    yip1_5th = yi + (2825*k1/27648 + 18575*k3/48384 + 13525*k4/55296 + 277*k5/14336 + k6/4)*h
    yerr = yip1_5th - yip1
    return yip1, yerr

# Chapra Fig. 25.22
def ode_adapt_1step(dfunc, xi, yi, htry, yscal):
    safety = 0.9
    econ = 1.89e-4
    EPS = 0.00005 # use machine precision?
    h = htry
    emax = np.nan
    while True:
        ytemp, yerr = ode_rk45_1step(dfunc, xi, yi, h)
        print("yerr = ", yerr)
        emax = abs(yerr/yscal/EPS)
        if emax <= 1:
            break
        # fix h
        htemp = safety*h*emax**(-0.25)
        h = max( abs(htemp), 0.25*abs(h) )
        print("Fix h = ", h)
    #
    if emax > econ:
        h_next = safety * emax**(-0.2) * h
    else:
        h_next = 4*h
    #
    x = xi + h
    y = ytemp
    return x, y, h_next


def ode_adapt(dfunc, xi, yi, htry, xf, NstepMax):
    SMALL = 1e-30
    #
    x = np.zeros(Nstep+1)
    y = np.zeros(Nstep+1)
    #
    x[0] = x0
    y[0] = y0
    h = htry
    NstepActual = Nstep
    for i in range(0,Nstep):
        print("\nstep = ", i + 1)
        #
        if x[i] + h > xf:
            h = xf - x[i]
        #
        print("Using h = ", h)
        dy = dfunc(x[i], y[i])
        yscal = abs(y[i]) + abs(h*dy) + SMALL
        print("yscal = ", yscal)
        x[i+1], y[i+1], h = ode_adapt_1step(dfunc, x[i], y[i], h, yscal)
        print("Solution: %18.10f %18.10f" % (x[i+1], y[i+1]))
        print("h_next = ", h)
        #
        if x[i+1] >= xf:
            NstepActual = i + 1
            break
    #
    return x[:NstepActual+1], y[:NstepActual+1]

def deriv(x, y):
    return 10*np.exp( -(x-2)**2/(2*0.075**2) ) - 0.6*y

def exact_sol(x):
    return 0.5*np.exp(-0.6*x)

x0 = 0.0
y0 = 0.5
xf = 4.0
h = 0.1
Nstep = int(xf/h)

import matplotlib.pyplot as plt
plt.clf()

x, y = ode_rk4(deriv, x0, y0, h, Nstep)
plt.plot(x, y, marker="o", label="RK4")

x, y = ode_adapt(deriv, x0, y0, 0.5, xf, Nstep)
plt.plot(x, y, marker="o", label="Adaptive")

plt.grid(True)
plt.savefig("IMG_chapra_example_25_14.pdf")


