from math import pi
from root_regula_falsi import *

r = 1.0
ρs = 200.0
ρw = 1000.0

def f(h):
    Vw = 4*pi*r**3/3 - pi*h**2/3*(3*r - h) # displaced volume of water
    Vs = 4*pi*r**3/3
    return ρw*Vw - ρs*Vs


xr = root_regula_falsi(f, 0.0, 2*r)