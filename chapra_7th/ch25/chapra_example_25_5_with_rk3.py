from math import exp

# 3rd order Runge-Kutta method (Eq. 25.39)
def ode_rk3_1step(dfunc, xi, yi, h):
    k1 = dfunc(xi, yi)
    k2 = dfunc(xi + 0.5*h, yi + 0.5*k1*h)
    k3 = dfunc(xi + h, yi - k1*h + 2*k2*h)
    yip1 = yi + (k1 + 4*k2 + k3)*h/6
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
    yp1 = ode_rk3_1step(deriv, x, y, h)
    y_true = exact_sol(xp1)
    ε_t = (y_true - yp1)/y_true * 100
    print("%f %12.7f %12.7f  %5.2f%%" % (xp1, y_true, yp1, abs(ε_t)))
    # For the next step
    x = xp1
    y = yp1
