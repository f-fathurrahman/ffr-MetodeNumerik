from sympy import *

# heat eq:
# u_t = α * u_xx + f
# -> f = u_t - α * u_xx

x, y, t = symbols("x y t")

r2 = x**2 + y**2
u = exp(-t)*exp(-r2)

nabla2u = simplify(diff(u, x, 2) + diff(u, y, 2))

f = simplify(diff(u, t) - nabla2u)
pprint(f)

