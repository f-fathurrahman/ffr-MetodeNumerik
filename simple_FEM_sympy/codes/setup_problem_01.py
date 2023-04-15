from sympy import *

init_printing(use_unicode=True)

# Spatial variable
x = Symbol("x", real=True)

# Domain length
L = Symbol("L", real=True, positive=True)


# Setup linear element
Nelements = 2 
Nnodes = Nelements + 1

# Setup the grid
xnodes = [Integer(0), L/3, L] # I use Integer(0) to make sure that all xnodes are SymPy objects
h = [xnodes[i+1]-xnodes[i] for i in range(Nnodes-1)]

assert(len(xnodes) == Nnodes)
assert(len(h) == Nelements)


# Setup basis functions
Nfuncs = []

# First
cond1 = (x >= xnodes[0]) & (x <= xnodes[1])
f1 = (xnodes[1] - x)/h[0]
Nfuncs.append( Piecewise( (f1, cond1), (0, True) ) )

# Second
cond1 = (x >= xnodes[0]) & (x <= xnodes[1])
f1 = ( x - xnodes[0])/h[0]
#
cond2 = (x >= xnodes[1]) & (x <= xnodes[2])
f2 = (xnodes[2] - x)/h[1]
#
Nfuncs.append( Piecewise( (f1, cond1), (f2, cond2), (0, True) ) )

# Third
cond1 = (x >= xnodes[Nnodes-2]) & (x <= xnodes[Nnodes-1])
f1 = (x - xnodes[Nnodes-2])/h[Nnodes-2]
Nfuncs.append( Piecewise( (f1, cond1), (0, True) ) )

