from sympy import *

init_printing(use_unicode=True)

T_L = Symbol("T_L", real=True, positive=True)
L = Symbol("L", real=True, positive=True)
q = Symbol("q", real=True, positive=True)
k = Symbol("k", real=True, positive=True)

x, y, z = symbols("x y z")

α = Symbol("alpha", real=True, positive=True)
A = Symbol("A", real=True, positive=True)
#Q = 1 # constant
Q = A*exp( -α*(z-L/2)**2 ) # make it as a function of z
pprint(Q)

expr1 = simplify( integrate(Q, (z,0,y)) )
pprint(expr1)

expr2 = simplify( integrate(expr1, (y,x,L)) )
pprint(expr2)

T_analytic = T_L + q/k*(L - x) + 1/k * expr2

q_L = simplify( k*diff(T_analytic,x).subs({x: L}) ) # positive value: q_L is pointing outside

from sympy.utilities.codegen import codegen

res = codegen( ("T_analytic", T_analytic), language="Julia")
print(res[0][1])

res = codegen( ("q_L", q_L), language="Julia")
print(res[0][1])
