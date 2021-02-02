def f(x):
    return -0.1*x**4 - 0.15*x**3 - 0.5*x**2 - 0.25*x + 1.2

def df(x):
    return -0.4*x**3 - 0.45*x**2 - x - 0.25

def centered_diff(f, x, h):
    return ( f(x+h) - f(x-h) )/(2*h)

x = 0.5
h = 1.0
true_val = df(x)

print("--------------------------------------------------------")
print("           h             approx_val             error")
print("--------------------------------------------------------")
for i in range(11):
    approx_val = centered_diff(f, x, h)
    εt = abs(approx_val - true_val)
    print("%18.10f %18.14f %18.13f" % (h, approx_val, εt))
    h = h/10

