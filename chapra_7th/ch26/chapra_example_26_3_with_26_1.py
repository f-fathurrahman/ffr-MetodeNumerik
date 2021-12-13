import numpy as np

# One-step application of Heun's method for ODE
# Using iterative steps to determine y0ip1
# non-self starting
# also return Ec (Eq. 26.26)
def ode_heun_1step_iterative(dfunc, xi, yi, yim1, h, NiterMax=100, Δ=1e-6):
    y0ip1 = yim1 + dfunc(xi,yi)*2*h # predictor
    #print("predictor: ", y0ip1)
    y0ip1_old = y0ip1
    xip1 = xi + h
    y0ip1_initial = y0ip1 # for calculating Ec
    for i in range(NiterMax+1):
        avg = 0.5*( dfunc(xi,yi) + dfunc(xip1,y0ip1) )*h
        y0ip1 = yi + avg
        diff = abs(y0ip1 - y0ip1_old)
        # Uncomment this to see the iteration process
        #print("iter: %2d y0ip1 = %12.7f  diff = %12.7e" % (i+1, y0ip1, diff))
        if diff <= Δ:
            break
        y0ip1_old = y0ip1
    Ec = -(y0ip1_initial - y0ip1)/5
    return y0ip1, Ec

def deriv(t, y):
    return -1000*y + 3000 - 2000*np.exp(-t)

# exact solution for y(0) = 0
def exact_sol(t):
    return 3 - 0.998*np.exp(-1000*t) - 2.002*np.exp(-t)


# For starting Heun's method
def ode_rk4_1step(dfunc, xi, yi, h):
    k1 = dfunc(xi, yi)
    k2 = dfunc(xi + 0.5*h, yi + 0.5*k1*h)
    k3 = dfunc(xi + 0.5*h, yi + 0.5*k2*h)
    k4 = dfunc(xi + h, yi + k3*h)
    yip1 = yi + (k1 + 2*k2 + 2*k3 + k4)*h/6
    return yip1


# initial cond
t0 = 0.0
y0 = 0.0

h = 0.0015
tf = 0.4
N = int((tf - t0)/h)
# Solution arrays
t = np.zeros(N+1)
y = np.zeros(N+1)

# Start with two points (intial + one point constructed from exact sol)
t[0] = t0
y[0] = y0

t[1] = t0 + h
#y[1] = exact_sol(t[1])
# alternative: using RK4
y[1] = ode_rk4_1step(deriv, t[0], y[0], h)

yim1 = y[0]
print("yim1 = ", yim1)
for i in range(0,N):
    t[i+1] = t[i] + h
    y[i+1], Ec = ode_heun_1step_iterative(deriv, t[i], y[i], yim1, h, NiterMax=10)
    y_true = exact_sol(t[i+1])
    ε_a = (y_true - y[i+1]) # error
    print("%f %12.7f %12.7f  %5.2f  %5.2f" % (t[i+1], y_true, y[i+1], abs(ε_a), Ec))
    # For the next step
    yim1 = y[i]

import matplotlib.pyplot as plt

plt.clf()
plt.plot(t, y, label="non-self starting Heun")
plt.plot(t, exact_sol(t), label="exact")

plt.xlim(0.0, 0.05)
plt.grid(True)
plt.legend()
plt.savefig("IMG_26_1_v2.pdf")