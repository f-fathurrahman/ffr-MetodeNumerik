import numpy as np
import matplotlib.pyplot as plt

from sympy import *

init_printing()

# x is symbolic object
# xgrid, h: array (numeric)
# L: numeric
def build_shape_functions( x, xgrid, h ):
    Nnodes = len(xgrid)
    
    # First shape functions
    cond1 = (x >= xgrid[0]) & (x <= xgrid[1])
    f1 = (xgrid[1] - x)/h[0]
    N_first = Piecewise( (f1, cond1), (0, True) )
    
    shape_functions = [ N_first ]
    
    for i in range(1,Nnodes-1):
        #
        cond1 = (x >= xgrid[i-1]) & (x <= xgrid[i])
        f1 = ( x - xgrid[i-1])/h[i-1]
        #
        cond2 = (x >= xgrid[i]) & (x <= xgrid[i+1])
        f2 = (xgrid[i+1] - x)/h[i]
        N_i = Piecewise( (f1, cond1), (f2, cond2), (0, True) )
        #
        shape_functions.append(N_i)
    
    cond1 = (x >= xgrid[Nnodes-2]) & (x <= xgrid[Nnodes-1])
    f1 = (x - xgrid[Nnodes-2])/h[Nnodes-2]
    N_last = Piecewise( (f1, cond1), (0, True) )
    shape_functions.append(N_last)
    
    return shape_functions

Nelements = 4
Nnodes = Nelements + 1

L = Symbol("L", real=True, positive=True)
k = Symbol("k", real=True, positive=True)
Q = Symbol("Q", real=True, positive=True)

xgrid = [ i*L/(Nnodes-1) for i in range(Nnodes) ]
h = [ (xgrid[i+1] - xgrid[i]) for i in range(Nnodes-1) ]
print(xgrid)
print(h)

x = Symbol("x", positive=True, real=True)
N = build_shape_functions( x, xgrid, h )

# Nodal values
#T_i = [ Symbol("T_1", real=True), Symbol("T_2", real=True),  Symbol("T_3", real=True) ]

T = [ Symbol("T_"+str(i+1), real=True) for i in range(Nnodes) ]
pprint(T)

ss = 0
for j in range(Nnodes):
    ss = ss + diff( N[j], x ) * T[j]

print("Coeffs:")
for i in range(Nnodes):
    expr1 = k * diff( N[i], x ) * ss
    expr2 = integrate(expr1, (x,0,L))
    pprint( expr2.expand().factor() )
