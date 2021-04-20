from sympy import *
x = symbols("x", real=True)
T = symbols("T", cls=Function)

my_eqn = T(x).diff(x,2) - 0.15*T(x)

sols = dsolve( my_eqn, T(x), ics = { T(0): 240, T(10): 150 } )
