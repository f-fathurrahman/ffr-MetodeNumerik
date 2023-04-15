# Linear element, arbitrary number of nodes

from sympy import *

init_printing(use_unicode=True)

# Spatial variable
x = Symbol("x", real=True)

# Domain length
#L = Symbol("L", real=True, positive=True)
L = 1.0

# Setup quadratic element
Nelements = 1
NnodesPerElement = 3
Nnodes = (NnodesPerElement-1)*Nelements + 1

#import numpy as np
#gblNum = np.zeros(NnodesPerElement, Nelements, dtype=np.int)


# Setup the grid (same spacing)
xnodes = [Integer(0), L/2, L]
h = [L]

assert(len(xnodes) == Nnodes)
assert(len(h) == Nelements)

# Setup basis functions
Nfuncs = []

Nfuncs.append( 1 - 3*x/h[0] + 2*(x/h[0])**2 )
Nfuncs.append( 4*x/h[0]*(1 - x/h[0]) )
Nfuncs.append( x/h[0]*(2*(x/h[0]) - 1) )
