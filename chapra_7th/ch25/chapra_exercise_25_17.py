from math import sqrt

def deriv(t, y):
    assert y >= 0.0
    k = 0.06
    return -k*sqrt(y)

# One-step application of Euler's method for ODE
def ode_euler_1step(dfunc, xi, yi, h):
    return yi + dfunc(xi,yi)*h

# RK5 method for ODE
def ode_rk5_1step(dfunc, xi, yi, h):
    k1 = dfunc(xi, yi)
    k2 = dfunc(xi + h/4, yi + k1*h/4)
    k3 = dfunc(xi + h/4, yi + k1*h/8 + k2*h/8)
    k4 = dfunc(xi + h/2, yi - k2*h/2 + k3*h)
    k5 = dfunc(xi + 3*h/4, yi + 3*k1*h/16 + 9*k4*h/16)
    k6 = dfunc(xi + h, yi - 3*k1*h/7 + 2*k2*h/7 + 12*k3*h/7 - 12*k4*h/7 + 8*k5*h/7)    
    yip1 = yi + (7*k1 + 32*k3 + 12*k4 + 32*k5 + 7*k6)*h/90
    return yip1

# Initial cond
t0 = 0.0
y0 = 3.0
yf = 0.0
h = 0.5 # try to change this

t = t0
y = y0
NstepMax = 1000
SMALL = 1e-5
print("%18.10f %18.10f" % (t, y))
#for i in range(0,NstepMax): # or use while True
while True:
    tp1 = t + h
    yp1 = ode_euler_1step(deriv, t, y, h)
    #yp1 = ode_rk5_1step(deriv, t, y, h)
    print("%18.10f %18.10f" % (tp1, yp1))
    if yp1 <= SMALL:
        break
    # For the next step
    t = tp1
    y = yp1
print("t = ", tp1)