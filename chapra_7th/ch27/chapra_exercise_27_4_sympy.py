from sympy import *
x = symbols("x", real=True)
y = symbols("y", cls=Function)

my_eqn = 7*y(x).diff(x,2) - 2*y(x).diff(x) - y(x) + x

sols = dsolve( my_eqn, y(x), ics = { y(0): 5, y(20): 8 } )
print(simplify(sols.args[1]))
pprint(sols)