from math import exp

# 5th order Runge-Kutta method (Eq. 25.41)
def ode_rk5_1step(dfunc, xi, yi, h):
    k1 = dfunc(xi, yi)
    k2 = dfunc(xi + h/4, yi + k1*h/4)
    k3 = dfunc(xi + h/4, yi + k1*h/8 + k2*h/8)
    k4 = dfunc(xi + h/2, yi - k2*h/2 + k3*h)
    k5 = dfunc(xi + 3*h/4, yi + 3*k1*h/16 + 9*k4*h/16)
    k6 = dfunc(xi + h, yi - 3*k1*h/7 + 2*k2*h/7 + 12*k3*h/7 - 12*k4*h/7 + 8*k5*h/7)    
    yip1 = yi + (7*k1 + 32*k3 + 12*k4 + 32*k5 + 7*k6)*h/90

    return yip1

def deriv(x, y):
    return 4*exp(0.8*x) - 0.5*y

def exact_sol(x):
    return 4.0/1.3*( exp(0.8*x) - exp(-0.5*x) ) + 2*exp(-0.5*x)

# Initial cond
x0 = 0.0
y0 = 2.0
xf = 4.0
h = 1.0
Nstep = int(xf/h)

x = x0
y = y0
for i in range(0,Nstep):
    xp1 = x + h
    yp1 = ode_rk5_1step(deriv, x, y, h)
    y_true = exact_sol(xp1)
    ε_t = (y_true - yp1)/y_true * 100
    print("%f %12.7f %12.7f  %8.4f%%" % (xp1, y_true, yp1, abs(ε_t)))
    # For the next step
    x = xp1
    y = yp1
