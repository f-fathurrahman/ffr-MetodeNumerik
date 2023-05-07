from sympy import *

x, y = symbols("x y")
u = exp(-(x**2 + y**2))*sin(2*x*y)
g = cos(2*y)

f = simplify(diff(u, x, 2) + diff(u, y, 2) + g*u)
pprint(f)

