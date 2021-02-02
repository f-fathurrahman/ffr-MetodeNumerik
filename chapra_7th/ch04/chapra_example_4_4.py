def f(x):
    return -0.1*x**4 - 0.15*x**3 - 0.5*x**2 - 0.25*x + 1.25

def df(x):
    return -0.4*x**3 - 0.45*x**2 - x - 0.25

x = 0.5

h = 0.5
#h = 0.25

true_val = df(x)
print("true_val = ", true_val)

print("h = ", h)

# Forward diff
approx_val = ( f(x+h) - f(x) )/h
print("Using forward diff: ", approx_val)
print("ε_t (in percent) = ", (true_val - approx_val)/true_val*100)

# Backward diff
approx_val = ( f(x) - f(x-h) )/h
print("Using backward diff: ", approx_val)
print("ε_t (in percent) = ", (true_val - approx_val)/true_val*100)

# Backward diff
approx_val = ( f(x+h) - f(x-h) )/(2*h)
print("Using centered diff: ", approx_val)
print("ε_t (in percent) = ", (true_val - approx_val)/true_val*100)
