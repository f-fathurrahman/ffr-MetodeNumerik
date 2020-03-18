from sympy import *

init_printing()

x = Symbol("x", real=True)
L = Symbol("L", real=True, positive=True)
k = Symbol("k", real=True, positive=True)
Q = Symbol("Q", real=True, positive=True)

T = Function("T")(x)
W = Function("W")(x)

diffEq = Equality( -k * Derivative(T,x,2), Q )
diffEq2 = Equality( W*diffEq.args[0], W*diffEq.args[1] )
#pprint(diffEq2)

LHS = -k * Derivative(T,x,2) - Q
pprint( W * LHS )