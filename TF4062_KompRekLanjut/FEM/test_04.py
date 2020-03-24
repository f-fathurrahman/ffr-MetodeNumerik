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

Nelements = 2
Nnodes = Nelements + 1
NnodesPerElement = 2 # linear element

L = Symbol("L", real=True, positive=True)
k = Symbol("k", real=True, positive=True)
Q = Symbol("Q", real=True, positive=True)

xgrid = [ i*L/(Nnodes-1) for i in range(Nnodes) ]
h = [ (xgrid[i+1] - xgrid[i]) for i in range(Nnodes-1) ]

x = Symbol("x", positive=True, real=True)
N = build_shape_functions( x, xgrid, h )

# Nodal values of
T_i = [ Symbol("T_1", real=True), Symbol("T_2", real=True),  Symbol("T_3", real=True) ]
expr1 = T_i[0] * N[0] + T_i[1] * N[1]

ss = 0
for j in range(Nnodes):
    ss = ss + diff(N[j], x) * T_i[j]

print("Coeffs:")
for i in range(NnodesPerElement):
    expr1 = diff(N[i],x) * ss * k
    expr2 = integrate( expr1, (x,0,L/2) )
    pprint( expr2.expand().factor() ) # Trying to get a suitable form

print("Q contrib")
for i in range(NnodesPerElement):
    expr1 = integrate( Q*N[i], (x,0,L/2) )
    pprint( expr1 )
