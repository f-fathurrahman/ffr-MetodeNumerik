import numpy as np

β_AdamsMoulton = {
    2: [1/2, 1/2],
    3: [5/12, 8/12, -1/12],
    4: [9/24, 19/24, -5/24, 1/24],
    5: [251/720, 646/720, -264/720, 106/720, -19/720],
    6: [475/1440, 1427/1440, -798/1440, 482/1440, -173/1440, 27/1440]
}


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


def ode_AdamsMoulton(dfunc, x0, y0, h, Nstep, Norder=3):
    β = β_AdamsMoulton[Norder]
    print(β)
    x = np.zeros(Nstep+1)
    y = np.zeros(Nstep+1)
    #
    x[0] = x0
    y[0] = y0
    for i in range(0,Norder):
        #print("Using RK4, i+1 = ", i+1)
        x[i+1] = x[i] + h
        y[i+1] = ode_rk4_1step(dfunc, x[i], y[i], h)
    #
    for i in range(Norder,Nstep):
        x[i+1] = x[i] + h
        print("Using AdamsMoulton, i+1 = ", i+1)
        yip1_old = y[i+1]
        for ii in range(1,10):
            s = 0.0
            for k in range(0,Norder):
                #print("access idx: ", i+1-k)
                s += β[k]*dfunc(x[i+1-k], y[i+1-k])
            y[i+1] = y[i] + h*s
            diff = abs(y[i+1] - yip1_old)
            print("iter = %d, diff = %e" % (ii, diff))
            if diff < 1e-8:
                break
            yip1_old = y[i+1]
    return x, y



def deriv(x, y):
    return 4*np.exp(0.8*x) - 0.5*y

def exact_sol(x):
    return 4.0/1.3*( np.exp(0.8*x) - np.exp(-0.5*x) ) + 2*np.exp(-0.5*x)


# Initial cond
x0 = 0.0
y0 = 2.0
xf = 4.0
h = 0.2
Nstep = int(xf/h)

x, y = ode_AdamsMoulton(deriv, x0, y0, h, Nstep, Norder=5)
_, yrk4 = ode_RK4(deriv, x0, y0, h, Nstep)

for i in range(0,Nstep+1):
    y_exact = exact_sol(x[i])
    Δ_AM = abs(y[i] - y_exact)
    Δ_RK4 = abs(yrk4[i] - y_exact)
    print("%18.5e %18.5e" % (Δ_AM, Δ_RK4))

import matplotlib.pyplot as plt

plt.clf()
plt.plot(x, y, label="AM5")
plt.plot(x, yrk4, label="RK4")
plt.grid(True)
plt.legend()
plt.savefig("IMG_test_AM_v1.pdf")

plt.clf()
plt.plot(x, np.abs(y - exact_sol(x)), marker="o", label="AM5")
plt.plot(x, np.abs(yrk4 - exact_sol(x)), marker="o", label="RK4")
plt.grid(True)
plt.legend()
plt.title("Abs error")
plt.savefig("IMG_test_AM_v2.pdf")
