from math import exp
from root_regula_falsi import *
from root_bisection import *
from root_ridder import *

g = 9.81
t = 4.0
v = 36.0
m = 82.0

def f(c):
    return v - g*m/c*(1 - exp(-c*t/m))

def calc_v(c):
    return g*m/c*(1 - exp(-c*t/m))

cl = 3.0
cu = 5.0
cr = root_bisection(f, cl, cu)
print("v = ", calc_v(cr), " (should be close to ", v, ")")

cl = 3.0
cu = 5.0
cr = root_regula_falsi(f, cl, cu)
print("v = ", calc_v(cr), " (should be close to ", v, ")")

cl = 3.0
cu = 5.0
cr = root_ridder(f, cl, cu)
print("v = ", calc_v(cr), " (should be close to ", v, ")")
