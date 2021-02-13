from root_bisection import *
from root_regula_falsi import *
from root_regula_falsi_mod import *
from root_ridder import *

def f(x):
    return x**4 - 8*x**3 - 35*x**2 + 450*x - 1001

xl = 5.5 #4.5
xu = 10.0 #6.0

xr = root_bisection(f, xl, xu)

xr = root_regula_falsi(f, xl, xu)

xr = root_regula_falsi_mod(f, xl, xu)

xr = root_ridder(f, xl, xu)
