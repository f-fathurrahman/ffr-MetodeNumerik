import numpy as np

def ode_rk4_1step(dfunc, xi, yi, h):
    k1 = dfunc(xi, yi)
    k2 = dfunc(xi + 0.5*h, yi + 0.5*k1*h)
    k3 = dfunc(xi + 0.5*h, yi + 0.5*k2*h)
    k4 = dfunc(xi + h, yi + k3*h)
    yip1 = yi + (k1 + 2*k2 + 2*k3 + k4)*h/6
    return yip1

def ode_rk4(dfunc, x0, y0, h, Nstep):
    x = np.zeros(Nstep+1)
    y = np.zeros(Nstep+1)
    #
    x[0] = x0
    y[0] = y0
    for i in range(0,Nstep):
        x[i+1] = x[i] + h
        y[i+1] = ode_rk4_1step(dfunc, x[i], y[i], h)
    return x, y

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
    return yip1, yerr

# Chapra Fig. 25.22
def ode_adapt_1step(dfunc, xi, yi, htry, yscal, EPS=1e-6):
    safety = 0.9
    econ = 1.89e-4
    h = htry
    emax = np.nan
    while True:
        ytemp, yerr = ode_rk45_1step(dfunc, xi, yi, h)
        emax = abs(yerr/yscal/EPS)
        if emax <= 1:
            break
        # fix h
        htemp = safety*h*emax**(-0.25)
        h = max( abs(htemp), 0.25*abs(h) )
    #
    if emax > econ:
        h_next = safety * emax**(-0.2) * h
    else:
        h_next = 4*h
    #
    x = xi + h
    y = ytemp
    return x, y, h_next


def ode_adapt(dfunc, xi, yi, htry, xf, NstepMax=100):
    SMALL = 1e-30
    #
    x = np.zeros(NstepMax+1)
    y = np.zeros(NstepMax+1)
    #
    x[0] = x0
    y[0] = y0
    h = htry
    NstepActual = 0
    for i in range(0,NstepMax):
        #
        if x[i] + h > xf:
            h = xf - x[i]
        #
        dy = dfunc(x[i], y[i])
        yscal = abs(y[i]) + abs(h*dy) + SMALL
        x[i+1], y[i+1], h = ode_adapt_1step(dfunc, x[i], y[i], h, yscal)
        #
        if x[i+1] >= xf:
            NstepActual = i + 1
            break
    #
    if x[NstepActual] < xf:
        print("WARNING: NstepMax is not sufficient x = ", x[NstepActual])
    return x[:NstepActual+1], y[:NstepActual+1]

def deriv(x, y):
    return 1/((x-0.3)**2 + 0.01) + 1/( (x - 0.9)**2 + 0.04 ) - 6

x0 = 0.0
y0 = 0.0
xf = 1.0
h = 0.01
for h in [0.5, 0.25, 0.1, 0.05, 0.01]:
    Nstep = int(xf/h)
    x, y = ode_rk4(deriv, x0, y0, h, Nstep)
    print("h = %18.10f Last points: %18.10f %18.10f" % (h, x[-1], y[-1]))

x_adapt, y_adapt = ode_adapt(deriv, x0, y0, 0.5, xf, NstepMax=100)
print("Adaptive Last points: %18.10f %18.10f" % (x[-1], y[-1]))


# Using Romberg integration

def integ_trapz_multiple( f, a, b, N ):
    x0 = a
    xN = b
    h = (b - a)/N
    ss = 0.0
    for i in range(1,N):
        xi = x0 + i*h
        ss = ss + f(xi)
    I = h/2 * ( f(x0) + 2*ss + f(xN) )
    return I

def integ_romberg(f, a, b, es=1e-10, MAXIT=10):
    I = np.zeros( (MAXIT+2,MAXIT+2) )
    n = 1
    # We start from I[1,1], to follow the book's notation
    I[1,1] = integ_trapz_multiple(f, a, b, n)
    iterConv = 0
    for i in range(1,MAXIT+1):
        n = 2**i
        I[i+1,1] = integ_trapz_multiple(f, a, b, n)
        #
        for k in range(2,i+2):
            j = 2 + i - k
            I[j,k] = ( 4**(k-1)*I[j+1,k-1] - I[j,k-1] )/ (4**(k-1) - 1)
        #
        ea = abs( (I[1,i+1] - I[2,i])/I[1,i+1] )*100 # in percent
        if ea <= es:
            iterConv = i
            print("converged, iterConv = ", iterConv)
            break
        iterConv = i
    print("iterConv = ", iterConv)
    # to make sure that we are use variable that is defined outside the loop
    # we use iterConv instead of i
    return I[1,iterConv+1]

def my_func(x):
    return 1/((x-0.3)**2 + 0.01) + 1/( (x - 0.9)**2 + 0.04 ) - 6

a = 0.0
b = 1.0
resN = integ_romberg(my_func, a, b)
print("integ_romberg = %18.10f" % resN)