from math import pi
from root_regula_falsi import *

R = 3.0  # m
V = 30.0 # m^3

def f(h):
    return V - pi*h**2*(3*R - h)/3

xroot = root_regula_falsi(f, 0.0, 3.0)
print(f(xroot))