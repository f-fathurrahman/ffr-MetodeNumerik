import numpy as np

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


def ode_Milne_1step(dfunc, xi, yi, yim1, yim2, yim3, h, NiterMax=100, Δ=1e-6):
    #print("yi   = ", yi)
    #print("yim1 = ", yim1)
    #print("yim2 = ", yim2)
    #print("yim3 = ", yim3)
    #
    fi = deriv(xi, yi)
    #
    xim1 = xi - h
    fim1 = deriv(xim1,yim1)
    #
    xim2 = xi - 2*h
    fim2 = deriv(xim2,yim2)
    #
    yip1 = yim3 + 4*h/3*( 2*fi - fim1 + 2*fim2 ) # predictor
    #print("fi = ", fi)
    #print("fim1 = ", fim1)
    #print("fim2 = ", fim2)
    #print("Predictor: ", yip1)
    yip1_old = yip1
    xip1 = xi + h
    for i in range(NiterMax+1):
        fip1 = deriv(xip1, yip1)
        yip1 = yim1 + h/3 * (fim1 + 4*fi + fip1)
        diff = abs(yip1 - yip1_old)
        # Uncomment this to see the iteration process
        #print("Milne iter: %2d yip1 = %12.7f  diff = %12.7e" % (i+1, yip1, diff))
        if diff <= Δ:
            break
        yip1_old = yip1
    return yip1



def deriv(x, y):
    return 4*np.exp(0.8*x) - 0.5*y

def exact_sol(x):
    return 4.0/1.3*( np.exp(0.8*x) - np.exp(-0.5*x) ) + 2*np.exp(-0.5*x)


# Initial cond
x0 = 0.0
y0 = 2.0
xf = 4.0
h = 1.0
Nstep = int(xf/h)

x = np.zeros(Nstep+1)
y = np.zeros(Nstep+1)

x[0] = x0
y[0] = y0

# Initial data, required for Milne algorithm
yim1 = exact_sol(x0-h)
yim2 = exact_sol(x0-2*h)
yim3 = exact_sol(x0-3*h)
x[1] = x[0] + h
y[1] = ode_Milne_1step(deriv, x[0], y[0], yim1, yim2, yim3, h)

print("%f %12.7f  %5.2f" % (x[1], y[1], abs(y[1]-exact_sol(x[1]))))

# Second step
yim3 = yim2
yim2 = yim1
yim1 = y[0]
x[2] = x[1] + h
y[2] = ode_Milne_1step(deriv, x[1], y[1], yim1, yim2, yim3, h)

print("%f %12.7f  %5.2f" % (x[2], y[2], abs(y[2]-exact_sol(x[2]))))


# Third step
yim3 = yim1
yim2 = y[0]
yim1 = y[1]
x[3] = x[2] + h
y[3] = ode_Milne_1step(deriv, x[2], y[2], yim1, yim2, yim3, h)

print("%f %12.7f  %5.2f" % (x[3], y[3], abs(y[3]-exact_sol(x[3]))))


# 4th step
yim3 = y[0]
yim2 = y[1]
yim1 = y[2]
x[4] = x[3] + h
y[4] = ode_Milne_1step(deriv, x[3], y[3], yim1, yim2, yim3, h)

print("%f %12.7f  %5.2f" % (x[4], y[4], abs(y[4]-exact_sol(x[4]))))
