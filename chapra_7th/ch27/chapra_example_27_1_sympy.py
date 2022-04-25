from sympy import *
x = symbols("x")
T = symbols("T", cls=Function)

Ta = 20
T1 = 40
T2 = 200
h = Rational(1,100) # 0.01
x1 = 0
x2 = 10

my_eqn = T(x).diff(x,2) + h*(Ta - T(x))

sols = dsolve( my_eqn, T(x), ics = { T(x1): T1, T(x2): T2 } )
print(sols)
pprint(sols)