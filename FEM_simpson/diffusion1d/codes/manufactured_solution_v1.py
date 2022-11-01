from sympy import *

x, t = symbols("x t", real=True)
H = symbols("L", real=True, positive=True)
κ = symbols("kappa")
L = symbols("L")
#u = x*(L - x)*5*t
u = x*(L - x)*exp(-t)

pde1 = diff(u, t) - κ*diff(u, x, 2)
pprint(simplify(pde1))