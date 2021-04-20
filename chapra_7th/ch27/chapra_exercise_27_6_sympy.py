from sympy import *
x = symbols("x", real=True)
T = symbols("T", cls=Function)

C1 = 1e-7
my_eqn = T(x).diff(x,2) - C1*T(x)**4 + 4*( 150.0 - (T(x) + 273.0) )

sols = dsolve( my_eqn, T(x), ics = { T(0): 200, T(0.5): 100 } )
print(simplify(sols.args[1]))
pprint(sols)