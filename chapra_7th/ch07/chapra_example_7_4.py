from scipy import optimize
import math

def f(x):
    return x - math.cos(x)

def fprime(x):
    return 1 + math.sin(x)

#sol = optimize.root_scalar(f, x0=0.0, fprime=fprime, method='newton')
#sol = optimize.root_scalar(f, x0=0.0, x1=1.0, method='secant')
#sol = optimize.root_scalar(f, bracket=[0.0, 1.0], method='brentq')
sol = optimize.root_scalar(f, bracket=[0.0, 1.0], method='bisect')
print(sol)
