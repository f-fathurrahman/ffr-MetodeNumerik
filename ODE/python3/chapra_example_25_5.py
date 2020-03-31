from math import exp

# One-step application of Heun's method for ODE
def ode_heun_1step(dfunc, xi, yi, h):
    y0ip1 = yi + dfunc(xi,yi)*h
    avg = 0.5*( dfunc(xi,yi) + dfunc(xi+h,y0ip1) )*h
    return yi + avg

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
    yp1 = ode_heun_1step(deriv, x, y, h)
    y_true = exact_sol(xp1)
    ε_t = (y_true - yp1)/y_true * 100
    print("%f %12.7f %12.7f  %5.2f%%" % (xp1, y_true, yp1, abs(ε_t)))
    # For the next step
    x = xp1
    y = yp1
