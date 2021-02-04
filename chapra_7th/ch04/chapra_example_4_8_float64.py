import numpy as np

def f(x):
    return -np.float64(0.1)*x**np.float64(4) - np.float64(0.15)*x**np.float64(3) - \
    np.float64(0.5)*x**np.float64(2) - np.float64(0.25)*x + np.float64(1.2)

def df(x):
    return -np.float64(0.4)*x**np.float64(3) - np.float64(0.45)*x**np.float64(2) - x - np.float64(0.25)

def centered_diff(f, x, h):
    return ( f(x+h) - f(x-h) )/(np.float64(2)*h)

x = np.float64(0.5)
h = np.float64(1.0)
true_val = df(x)

print("--------------------------------------------------------")
print("           h             approx_val             error")
print("--------------------------------------------------------")

for i in range(11):
    approx_val = centered_diff(f, x, h)
    εt = abs(approx_val - true_val)
    print("%18.10f %18.14f %18.13f" % (h, approx_val, εt))
    h = h/np.float64(10)

print(type(h))
print(type(centered_diff(f,x,h)))
