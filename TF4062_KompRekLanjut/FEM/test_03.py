import numpy as np
import matplotlib.pyplot as plt

from sympy import *

init_printing()

# x is symbolic object
# xgrid, h: array (numeric)
# L: numeric
def build_shape_functions( x, xgrid, h ):
    Npoints = len(xgrid)
    
    # First shape functions
    cond1 = (x >= xgrid[0]) & (x <= xgrid[1])
    f1 = (xgrid[1] - x)/h[0]
    N_first = Piecewise( (f1, cond1), (0, True) )
    
    shape_functions = [ N_first ]
    
    for i in range(1,Npoints-1):
        #
        cond1 = (x >= xgrid[i-1]) & (x <= xgrid[i])
        f1 = ( x - xgrid[i-1])/h[i-1]
        #
        cond2 = (x >= xgrid[i]) & (x <= xgrid[i+1])
        f2 = (xgrid[i+1] - x)/h[i]
        N_i = Piecewise( (f1, cond1), (f2, cond2), (0, True) )
        #
        shape_functions.append(N_i)
    
    cond1 = (x >= xgrid[Npoints-2]) & (x <= xgrid[Npoints-1])
    f1 = (x - xgrid[Npoints-2])/h[Npoints-2]
    N_last = Piecewise( (f1, cond1), (0, True) )
    shape_functions.append(N_last)
    
    return shape_functions

Nelements = 2
Npoints = Nelements + 1

L = Symbol("L", real=True, positive=True)
k = Symbol("k", real=True, positive=True)
Q = Symbol("Q", real=True, positive=True)

#L = 1
xgrid = [ i*L/(Npoints-1) for i in range(Npoints) ] #[0, L/2, L]
h = [ (xgrid[i+1] - xgrid[i]) for i in range(Npoints-1) ]
print(xgrid)
print(h)

x = Symbol("x", positive=True, real=True)
N = build_shape_functions( x, xgrid, h )

T = Function("T")(x)
W = Function("W")(x)

LHS = -k * Derivative(T,x,2) - Q

T_i = [ Symbol("T_1", real=True), Symbol("T_2", real=True),  Symbol("T_3", real=True) ]
expr1 = T_i[0] * N[0] + T_i[1] * N[1]

#pprint( LHS.subs({T: expr1}) )
#pprint( LHS*N[0] )

#pprint( LHS.subs)

ss = 0
for j in range(Npoints):
    ss = ss + diff(N[j], x) * T_i[j]

print("Coeffs:")
for i in range(Npoints):
    expr1 = diff(N[i],x) * ss * k
    #pprint(expr1)
    expr2 = integrate(expr1, (x,0,L/2)) + integrate(expr1, (x,L/2,L))
    pprint( expr2.expand().factor() )

pprint( integrate( N[0]*Q, (x,0,L/2) ) )