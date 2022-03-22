import numpy as np
from natural_cubic_spline import *

Npts_known = 7
x = np.zeros(Npts_known)
y = np.zeros(Npts_known)

L = 1.0
dx = L/(Npts_known-1)
print("dx = ", dx)
for i in range(0,Npts_known):
    x[i] = i*dx
    y[i] = np.cos(2.0*np.pi*x[i]/L)

#N = len(x) - 1
#allocate( e(N-1), f(N-1), g(N-1), r(N-1) )
#allocate( d2x(0:N) )
#! Natural spline
#d2x(0) = 0.d0
#d2x(N) = 0.d0

for i in range(0,Npts_known):
    print("%18.10f %18.10f" % (x[i], y[i]))

e, f, g, r = gen_trid_matrix( x, y)
print("After gen_trid_matrix")
print("e = ", e)
print("f = ", f)
print("g = ", g)
print("r = ", r)

decomp_trid(e, f, g)
print("After decomp_trid: ")
print("e = ", e)
print("f = ", f)
print("g = ", g)
print("r = ", r)

d2x = np.zeros(Npts_known)
# Natural spline
d2x[0] = 0.0
d2x[Npts_known-1] = 0.0 # last element
d2x[1:Npts_known-1] = subs_trid( e, f, g, r )

print("After subs_trid: ")
print("e = ", e)
print("f = ", f)
print("g = ", g)
print("r = ", r)
print("d2x = ", d2x)

xu = 0.51
y, dy, d2y = interp_nat_cubic_spline( x, y, d2x, xu )
print("yu  = ", y)
print("dy  = ", dy)
print("d2y = ", d2y)
