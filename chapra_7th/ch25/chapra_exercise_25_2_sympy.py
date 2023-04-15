# FIXME: Broken in SymPy v1.9 and v.10
# works with v1.7

from sympy import *
t = symbols("t", real=True)
y = symbols("y", cls=Function)

my_eqn = y(t).diff(t,1) - (1 + 4*t)*sqrt(y(t))

sols = dsolve( my_eqn, y(t) )

#sols = dsolve( my_eqn, y(t), ics = { y(0): 1 } )
print(sols)
pprint(sols)
