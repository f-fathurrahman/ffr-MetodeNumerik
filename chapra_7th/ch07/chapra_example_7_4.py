from scipy import optimize
import math

def f(x):
    return x - math.cos(x)

def fprime(x):
    return 1 + math.sin(x)

print("\nUsing Newton-Raphson method")
sol = optimize.root_scalar(f, x0=0.0, fprime=fprime, method='newton')
print(sol)

print("\nUsing secant method")
sol = optimize.root_scalar(f, x0=0.0, x1=1.0, method='secant')
print(sol)

# Bracketing method
for method in ["brentq", "brenth", "ridder", "bisect"]:
    print("\nUsing %s method" % (method))
    sol = optimize.root_scalar(f, bracket=[0.0, 1.0], method=method)
    print(sol)

print("\nUsing bisect directly")
xroot, sol = optimize.bisect(f, a=0.0, b=1.0, full_output=True)
print(sol)
