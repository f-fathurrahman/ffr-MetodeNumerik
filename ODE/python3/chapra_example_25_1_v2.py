# the left hand side of dy/dx=... (in general depends on x and y)
# In the present case it only depends on x
def deriv(x, y):
    return -2*x**3 + 12*x**2 - 20*x + 8.5

# One-step application of Euler's method for ODE
def ode_euler_1step(dfunc, xi, yi, h):
    return yi + dfunc(xi,yi)*h

def exact_sol(x):
    return -0.5*x**4 + 4*x**3 - 10*x**2 + 8.5*x + 1

# initial cond
x0 = 0.0
y0 = 1.0

print("   x     y_true   y_Euler       ε_t")
print("-------------------------------------")

print("%5.1f  %8.5f  %8.5f" % (x0, y0, y0)) # Intial cond

x = x0
y = y0
h = 0.5
for i in range(0,8):
    xp1 = x + h
    yp1 = ode_euler_1step(deriv, x, y, h)
    y_true = exact_sol(xp1)
    ε_t = (y_true - yp1)/y_true * 100
    print("%5.1f  %8.5f  %8.5f  %8.1f %%" % (xp1, y_true, yp1, ε_t))
    # Update x and y for the next step
    x = xp1
    y = yp1
