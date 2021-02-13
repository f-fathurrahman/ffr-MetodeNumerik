from math import exp
from root_regula_falsi import *
from root_bisection import *
from root_ridder import *

g = 9.81
c = 15
t = 10.0
v = 36.0

def f(m):
    return v - g*m/c*(1 - exp(-c*t/m))

def calc_v(m):
    return g*m/c*(1 - exp(-c*t/m))

ml = 1.0
mu = 100.0
mr = root_bisection(f, ml, mu)
print("v = ", calc_v(mr), " (should be close to ", v, ")")

ml = 1.0
mu = 100.0
mr = root_regula_falsi(f, ml, mu)
print("v = ", calc_v(mr), " (should be close to ", v, ")")

ml = 1.0
mu = 100.0
mr = root_ridder(f, ml, mu)
print("v = ", calc_v(mr), " (should be close to ", v, ")")
