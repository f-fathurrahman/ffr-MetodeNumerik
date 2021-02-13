from math import pi
from root_regula_falsi import *
from root_bisection import *
from root_regula_falsi_mod import *
from root_ridder import *

r1 = 0.5
r2 = 1.0
h = 1.0
ρf = 200.0
ρw = 1000.0


def calc_V_frustum(r1, r2, h):
    return pi*h/3*(r1**2 + r2**2 + r1*r2)

# h1/h = (r2-r1)/r
# -> r = (r2-r1)*h/h1

def f(h1):
    r = (r2 - r1)*h/h1
    Vw = calc_V_frustum(r1 + r, r2, h-h1) # displaced volume of water
    Vf = calc_V_frustum(r1, r2, h) # cone volume
    return ρw*Vw - ρf*Vf

xr = root_regula_falsi(f, 0.1, 1.0, NiterMax=1000)
xr = root_regula_falsi_mod(f, 0.1, 1.0, NiterMax=1000)
xr = root_bisection(f, 0.01, 1.0, NiterMax=1000)
xr = root_ridder(f, 0.01, 1.0, NiterMax=1000)