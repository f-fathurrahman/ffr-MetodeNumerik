import numpy as np

from natural_cubic_spline import *

x = np.array([3.0, 4.5, 7.0, 9.0])
y = np.array([2.5, 1.0, 2.5, 0.5])

e, f, g, r = gen_trid_matrix(x, y)
print("e = ", e) # subdiagonal
print("f = ", f) # main diagonal
print("g = ", g) # subdiagonal
print("r = ", r) # rhs vector

Npts_known = len(x)
d2x = np.zeros(len(x))
# Natural spline conditions
d2x[0] = 0.0
d2x[Npts_known-1] = 0.0 # last element

# Solve tridiagonal system (decomposition and back substitution)
decomp_trid(e, f, g)
d2x[1:Npts_known-1] = subs_trid( e, f, g, r )
print("d2x = ", d2x)

xu = 5
y, dy, d2y = interp_nat_cubic_spline( x, y, d2x, xu )
print("y = ", y)

