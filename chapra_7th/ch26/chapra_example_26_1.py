import numpy as np
import matplotlib.pyplot as plt

def deriv(t, y):
    return -1000*y + 3000 - 2000*np.exp(-t)

# explicit Euler
def ode_euler_1step(dfunc, ti, yi, h):
    return yi + dfunc(ti,yi)*h

# implicit Euler
# dy/dt is hardcoded
def ode_euler_1step_implicit(tip1, yi, h):
    num = yi + 3000*h - 2000*h*np.exp(-tip1)
    denum = 1 + 1000*h
    return num/denum

# exact solution for y(0) = 0
def exact_sol(t):
    return 3 - 0.998*np.exp(-1000*t) - 2.002*np.exp(-t)


plt.clf()

# initial cond
t0 = 0.0
y0 = 0.0

h = 0.0005
tf = 0.4
N = int((tf - t0)/h)
# Solution arrays
t = np.zeros(N+1)
y = np.zeros(N+1)
t[0] = t0
y[0] = y0
for i in range(0,N):
    t[i+1] = t[i] + h
    y[i+1] = ode_euler_1step(deriv, t[i], y[i], h)
plt.plot(t, y, label="explicit Euler h="+str(h))

# Plot exact solution, using the densest t grid, i.e. the current
plt.plot(t, exact_sol(t), label="exact")

h = 0.0015
tf = 0.4
N = int((tf - t0)/h)
# Solution arrays
t = np.zeros(N+1)
y = np.zeros(N+1)
t[0] = t0
y[0] = y0
for i in range(0,N):
    t[i+1] = t[i] + h
    y[i+1] = ode_euler_1step(deriv, t[i], y[i], h)
plt.plot(t, y, label="explicit Euler h=" + str(h))



h = 0.01
tf = 0.4
N = int((tf - t0)/h)
# Solution arrays
t = np.zeros(N+1)
y = np.zeros(N+1)
t[0] = t0
y[0] = y0
for i in range(0,N):
    t[i+1] = t[i] + h
    y[i+1] = ode_euler_1step_implicit(t[i+1], y[i], h)
plt.plot(t, y, marker="o", label="implicit Euler h=" + str(h))


plt.xlim(0.0, 0.05)
plt.grid(True)
plt.legend()
plt.savefig("IMG_chapra_example_36_1.pdf")
