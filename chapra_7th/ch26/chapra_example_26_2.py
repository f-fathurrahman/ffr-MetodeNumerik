from math import exp

# One-step application of Heun's method for ODE
# Using iterative steps to determine y0ip1
# non-self starting
def ode_heun_1step_iterative(dfunc, xi, yi, yim1, h, NiterMax=100, Δ=1e-6):
    y0ip1 = yim1 + dfunc(xi,yi)*2*h # predictor
    print("predictor: y0ip1 = ", y0ip1)
    y0ip1_old = y0ip1
    xip1 = xi + h
    for i in range(NiterMax+1):
        avg = 0.5*( dfunc(xi,yi) + dfunc(xip1,y0ip1) )*h
        y0ip1 = yi + avg
        diff = abs(y0ip1 - y0ip1_old)
        # Uncomment this to see the iteration process
        print("iter: %2d y0ip1 = %12.7f  diff = %12.7e" % (i+1, y0ip1, diff))
        if diff <= Δ:
            break
        y0ip1_old = y0ip1
    return y0ip1

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
#yim1 = -0.3929953 # value of solution at y(x=-1)
yim1 = exact_sol(-1.0)
for i in range(0,Nstep):
    xp1 = x + h
    yp1 = ode_heun_1step_iterative(deriv, x, y, yim1, h)
    y_true = exact_sol(xp1)
    ε_t = (y_true - yp1)/y_true * 100
    print("%f %12.7f %12.7f  %5.2f%%" % (xp1, y_true, yp1, abs(ε_t)))
    # For the next step
    yim1 = y
    x = xp1
    y = yp1
    print("x = ", x, " y = ", y, " yim1 = ", yim1)
    print()
