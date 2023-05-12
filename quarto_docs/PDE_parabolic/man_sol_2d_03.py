from sympy import *

init_printing()

# heat eq:
# u_t = α * u_xx + f
# -> f = u_t - α * u_xx

x, y, t = symbols("x y t")
Lx, Ly = symbols("L_x L_y")

u = 5*t*x*(Lx - x)*y*(y - Ly)

nabla2u = simplify(diff(u, x, 2) + diff(u, y, 2))

f = simplify(diff(u, t) - nabla2u)
pprint(f)

