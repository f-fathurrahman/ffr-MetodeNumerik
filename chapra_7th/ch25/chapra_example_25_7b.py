from math import exp

def deriv(x, y):
    return 4*exp(0.8*x) - 0.5*y

def exact_sol(x):
    return 4.0/1.3*( exp(0.8*x) - exp(-0.5*x) ) + 2*exp(-0.5*x)

# Classical 4th order Runge-Kutta method (Eq. 25.40)
def ode_rk4_1step(dfunc, xi, yi, h):
    k1 = dfunc(xi, yi)
    k2 = dfunc(xi + 0.5*h, yi + 0.5*k1*h)
    k3 = dfunc(xi + 0.5*h, yi + 0.5*k2*h)
    k4 = dfunc(xi + h, yi + k3*h)
    yip1 = yi + (k1 + 2*k2 + 2*k3 + k4)*h/6
    return yip1

# initial cond
x0 = 0.0
y0 = 2.0

# Using step of 0.5, starting from x0 and y0
x = x0
y = y0
h = 0.5
xp1 = x + h # we are searching for solution at x = 0.5
yp1 = ode_rk4_1step(deriv, x, y, h)
y_true = exact_sol(xp1)
ε_t = (y_true - yp1)/y_true * 100
print("First step : x = %f y_true = %.5f y = %.5f ε_t = %.4f %%" % (xp1, y_true, yp1, ε_t))
