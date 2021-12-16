from math import exp

def ode_rk4_1step(dfunc, xi, yi, h):
    k1 = dfunc(xi, yi)
    k2 = dfunc(xi + 0.5*h, yi + 0.5*k1*h)
    k3 = dfunc(xi + 0.5*h, yi + 0.5*k2*h)
    k4 = dfunc(xi + h, yi + k3*h)
    yip1 = yi + (k1 + 2*k2 + 2*k3 + k4)*h/6
    return yip1

def deriv(x, y):
    return 4*exp(0.8*x) - 0.5*y

def exact_sol(x):
    return 4.0/1.3*( exp(0.8*x) - exp(-0.5*x) ) + 2*exp(-0.5*x)

x0 = 0.0
xf = 2.0
h = 2.0
y0 = 2.0

y_exact = exact_sol(xf)

# single prediction with step h
y_h = ode_rk4_1step(deriv, x0, y0, h)

# two predictions, each with step h/2
y_h2 = ode_rk4_1step(deriv, x0, y0, h/2)
y_h2 = ode_rk4_1step(deriv, x0+h/2, y_h2, h/2)

print("y_h   = %18.10f error = %10.5e" % (y_h, abs(y_h - y_exact)))
print("y_h2  = %18.10f error = %10.5e" % (y_h2, abs(y_h2 - y_exact)))

Δ = y_h2 - y_h
ynew = y_h2 + Δ/15
print("y_new = %18.10f error = %10.5e" % (ynew, abs(ynew - y_exact)))

