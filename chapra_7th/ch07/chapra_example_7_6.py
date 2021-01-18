from scipy import optimize
from scipy.optimize import fsolve
import math

def f(x):
    return x**10 - 1

#sol = optimize.root_scalar(f, bracket=[-1.3, 0.0], method='bisect')
#sol = optimize.root_scalar(f, bracket=[-1.3, 0.0], method='ridder')

sol = optimize.root_scalar(f, bracket=[-1.3, 0.0], method='brentq')
print(sol)
sol = optimize.root_scalar(f, bracket=[0.0, 1.3], method='brentq')
print(sol)

