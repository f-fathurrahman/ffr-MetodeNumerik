# initial cond
x0 = 0.0
y0 = 1.0

# Using step of 0.5, starting from x0 and y0
x = x0
y = y0
h = 0.5
xp1 = x + h # we are searching for solution at x = 0.5
ϕ = -2*x**3 + 12*x**2 - 20*x + 8.5
yp1 = y + ϕ*h
y_true = -0.5*xp1**4 + 4*xp1**3 - 10*xp1**2 + 8.5*xp1 + 1
ε_t = (y_true - yp1)/y_true * 100
print("First step : x = %f y_true = %.5f y = %.5f ε_t = %.1f %%" % (xp1, y_true, yp1, ε_t))

# Second step
x = xp1 # x from the previous step
y = yp1 # y from the previoud step
h = 0.5
xp1 = x + h # we are searching for solution at x = 1.0
ϕ = -2*x**3 + 12*x**2 - 20*x + 8.5
yp1 = y + ϕ*h
y_true = -0.5*xp1**4 + 4*xp1**3 - 10*xp1**2 + 8.5*xp1 + 1
ε_t = (y_true - yp1)/y_true * 100
print("Second step: x = %f y_true = %.5f y = %.5f ε_t = %.1f %%" % (xp1, y_true, yp1, ε_t))