# FIXME: Broken in SymPy v1.9 and v.10
# works with v1.7

from sympy import *
x = symbols("x", real=True)
T = symbols("T", cls=Function)

my_eqn = T(x).diff(x,2) - 0.15*T(x)

sols = dsolve( my_eqn, T(x), ics = { T(0): 240, T(10): 150 } )
print(sols)
pprint(sols)