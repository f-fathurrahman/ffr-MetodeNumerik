import numpy as np

β_AdamsBashforth = {
    1: [1.0],
    2: [3/2, -1/2],
    3: [23/12, -16/12, 5/12],
    4: [55/24, -59/24, 37/24, -9/24],
    5: [1901/720, -2774/720, 2616/720, -1274/720, 251/720],
    6: [4277/720, -7923/720, 9982/720, -7298/720, 2877/720, -475/720]
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


def ode_AdamsBashforth(dfunc, x0, y0, h, Nstep, Norder=3):
    β = β_AdamsBashforth[Norder]
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
        s = 0.0
        #print("Using AdamsBashforth, i+1 = ", i+1)
        for k in range(0,Norder):
            print("i = ", i, "k = ", k)
            s += β[k]*dfunc(x[i-k], y[i-k])
        y[i+1] = y[i] + h*s

    return x, y



def deriv(x, y):
    return 4*np.exp(0.8*x) - 0.5*y

def exact_sol(x):
    return 4.0/1.3*( np.exp(0.8*x) - np.exp(-0.5*x) ) + 2*np.exp(-0.5*x)


# Initial cond
x0 = 0.0
y0 = 2.0
xf = 4.0
h = 0.1
Nstep = int(xf/h)

x, y = ode_AdamsBashforth(deriv, x0, y0, h, Nstep, Norder=2)
_, yrk4 = ode_RK4(deriv, x0, y0, h, Nstep)


for i in range(0,Nstep+1):
    y_exact = exact_sol(x[i])
    Δ_AB = abs(y[i] - y_exact)
    Δ_RK4 = abs(yrk4[i] - y_exact)
    print("%18.5e %18.5e %d" % (Δ_AB, Δ_RK4, Δ_AB < Δ_RK4) )
