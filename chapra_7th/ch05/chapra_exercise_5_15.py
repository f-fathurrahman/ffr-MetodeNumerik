from math import sqrt, tanh
from root_regula_falsi import *
from root_bisection import *
from root_ridder import *

g = 9.81
v = 5.0
t = 2.5
L = 4.0

def f(H):
    x = sqrt(2*g*H)
    return v - x*tanh(x/(2*L)*t)

def calc_v(H):
    x = sqrt(2*g*H)
    return x*tanh(x/(2*L)*t)

Hl = 0.0
Hu = 2.0
Hr = root_bisection(f, Hl, Hu)
print("v = ", calc_v(Hr), " (should be close to ", v, ")")

Hr = root_regula_falsi(f, Hl, Hu)
print("v = ", calc_v(Hr), " (should be close to ", v, ")")

Hr = root_ridder(f, Hl, Hu)
print("v = ", calc_v(Hr), " (should be close to ", v, ")")