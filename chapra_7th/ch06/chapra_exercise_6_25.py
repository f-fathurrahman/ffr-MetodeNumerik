from math import sqrt
from root_secant import *

scale = 1 #1e-5
W = 1e6*scale
k = 0.25
V = 1e6*scale
Q = 1e5*scale

def f(c):
    return W - Q*c - k*V*sqrt(c)

c0 = 4.0
δ = 0.5

# Need to use lower than usual to converge, may be due to round off error
# use scale of about 1e-5 for better convergence
cr = root_secant_mod(f, c0, δ=δ, TOL=1e-9)

print("try again with default δ")
cr = root_secant_mod(f, c0, TOL=1e-9)