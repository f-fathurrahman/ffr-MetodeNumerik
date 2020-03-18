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

Nelements = 3
Npoints = Nelements + 1

L = 1
xgrid = [ i*L/(Npoints-1) for i in range(Npoints) ] #[0, L/2, L]
h = [ (xgrid[i+1] - xgrid[i]) for i in range(Npoints-1) ]
print(xgrid)
print(h)

x = symbols("x")
N = build_shape_functions( x, xgrid, h )

NptsPlot = 100
xplt = np.linspace(0.0, L, NptsPlot)
yplt = np.zeros(NptsPlot)

plt.clf()

for i in range(Npoints):
    for ip in range(NptsPlot):
        yplt[ip] = N[i].subs( {x: xplt[ip]} )
    plt.plot(xplt, yplt, label="N_"+str(i+1))
plt.legend()
plt.savefig("IMG_ShapeFunctions.pdf")

#print(N[0].subs({x: 0.0}))

# Using SymPy's plot
#plot(N[0], N[1], N[2], (x,0,L))
#plot(N[2], (x,0,L))

