import numpy as np

def ode_rk4_1step(dfunc, xi, yi, h):
    k1 = dfunc(xi, yi)
    k2 = dfunc(xi + 0.5*h, yi + 0.5*k1*h)
    k3 = dfunc(xi + 0.5*h, yi + 0.5*k2*h)
    k4 = dfunc(xi + h, yi + k3*h)
    yip1 = yi + (k1 + 2*k2 + 2*k3 + k4)*h/6
    return yip1

def ode_RK4(dfunc, x0, y0, h, Nstep):
    x = np.zeros(Nstep+1)
    y = np.zeros(Nstep+1)
    #
    x[0] = x0
    y[0] = y0
    for i in range(0,Nstep):
        x[i+1] = x[i] + h
        y[i+1] = ode_rk4_1step(dfunc, x[i], y[i], h)
    return x, y


def ode_Adams4th_1step(dfunc, xi, yi, yim1, yim2, yim3, h, NiterMax=100, Δ=1e-6):
    fi = deriv(xi, yi)
    #
    xim1 = xi - h
    fim1 = deriv(xim1,yim1)
    #
    xim2 = xi - 2*h
    fim2 = deriv(xim2,yim2)
    #
    xim3 = xi - 3*h
    fim3 = deriv(xim3,yim3)
    #
    yip1 = yi + h/24 * ( 55*fi - 59*fim1 + 37*fim2 - 9*fim3 ) # predictor
    yip1_old = yip1
    xip1 = xi + h
    for i in range(NiterMax+1):
        fip1 = deriv(xip1, yip1)
        yip1 = yi + h/24 * (9*fip1 + 19*fi - 5*fim1 + fim2)
        diff = abs(yip1 - yip1_old)
        # Uncomment this to see the iteration process
        #print("Adams4th iter: %2d yip1 = %12.7f  diff = %12.7e" % (i+1, yip1, diff))
        if diff <= Δ:
            break
        yip1_old = yip1
    return yip1

def ode_Adams4th(dfunc, x0, y0, h, Nstep):
    #
    assert Nstep > 4
    #
    x = np.zeros(Nstep+1)
    y = np.zeros(Nstep+1)
    #
    x[0] = x0
    y[0] = y0
    #
    # Initial data, using ode_RK4
    for i in range(0,4):
        x[i+1] = x[i] + h
        y[i+1] = ode_rk4_1step(dfunc, x[i], y[i], h)
    #
    for i in range(4,Nstep):
        x[i+1] = x[i] + h
        y[i+1] = ode_Adams4th_1step(dfunc, x[i], y[i], y[i-1], y[i-2], y[i-3], h)
    return x, y


def ode_Milne_1step(dfunc, xi, yi, yim1, yim2, yim3, h, NiterMax=100, Δ=1e-6):
    #
    fi = deriv(xi, yi)
    #
    xim1 = xi - h
    fim1 = deriv(xim1,yim1)
    #
    xim2 = xi - 2*h
    fim2 = deriv(xim2,yim2)
    #
    yip1 = yim3 + 4*h/3*( 2*fi - fim1 + 2*fim2 ) # predictor
    yip1_old = yip1
    xip1 = xi + h
    for i in range(NiterMax+1):
        fip1 = deriv(xip1, yip1)
        yip1 = yim1 + h/3 * (fim1 + 4*fi + fip1)
        diff = abs(yip1 - yip1_old)
        # Uncomment this to see the iteration process
        #print("Milne iter: %2d yip1 = %12.7f  diff = %12.7e" % (i+1, yip1, diff))
        if diff <= Δ:
            break
        yip1_old = yip1
    return yip1


def ode_Milne(dfunc, x0, y0, h, Nstep):
    #
    assert Nstep > 4
    #
    x = np.zeros(Nstep+1)
    y = np.zeros(Nstep+1)
    #
    x[0] = x0
    y[0] = y0
    #
    # Initial data, using ode_RK4
    for i in range(0,4):
        x[i+1] = x[i] + h
        y[i+1] = ode_rk4_1step(dfunc, x[i], y[i], h)
    #
    for i in range(4,Nstep):
        x[i+1] = x[i] + h
        y[i+1] = ode_Milne_1step(dfunc, x[i], y[i], y[i-1], y[i-2], y[i-3], h)
    return x, y


def deriv(x, y):
    return -y

def exact_sol(x):
    return np.exp(-x)

# Initial cond
x0 = 0.0
y0 = 1.0
xf = 10.0
h = 0.5
Nstep = int(xf/h)

x, yAdams4th = ode_Adams4th(deriv, x0, y0, h, Nstep)
x, yMilne = ode_Milne(deriv, x0, y0, h, Nstep)

import matplotlib.pyplot as plt
plt.clf()
plt.plot(x, yAdams4th, marker="o", label="Adams4th")
plt.plot(x, yMilne, marker="o", label="Milne")
plt.plot(x, exact_sol(x), label="exact")
plt.grid(True)
plt.legend()
plt.xlim(5, 10)
plt.ylim(np.min(yMilne), 0.005)
plt.savefig("IMG_chapra_example_26_7.pdf")

for i in range(0,Nstep):
    Δ_Adams4th = np.abs(exact_sol(x[i]) - yAdams4th[i])
    Δ_Milne = np.abs(exact_sol(x[i]) - yMilne[i])
    print("%18.10e %18.10e" % (Δ_Adams4th, Δ_Milne))
