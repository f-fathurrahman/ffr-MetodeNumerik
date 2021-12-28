from math import exp

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
    print("k1 = %18.6f" % k1)
    print("k2 = %18.6f" % k2)
    print("k3 = %18.6f" % k3)
    print("k4 = %18.6f" % k4)
    print("k5 = %18.6f" % k5)
    print("k6 = %18.6f" % k6)
    print("yip1 = %18.6f" % yip1)
    print("yip1_5th = %18.6f" % yip1_5th)
    return yip1, yerr


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
y_h, yerr = ode_rk45_1step(deriv, x0, y0, h)
print(y_h)
print(yerr)

