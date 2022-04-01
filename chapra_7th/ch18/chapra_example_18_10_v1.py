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
print("d2x = ", d2x) # compare with Chapra's

xu = 5
# We do not use the derivatives, we use underscores to represent them
yu, _, _ = interp_nat_cubic_spline( x, y, d2x, xu )
print("yu = ", yu) # compare with Chapra's

# Make a plot
NptsPlot = 50
xgrid = np.linspace(3.0, 9.0, NptsPlot)
ygrid = np.zeros(NptsPlot)
for i in range(NptsPlot):
    ygrid[i], _, _ = interp_nat_cubic_spline( x, y, d2x, xgrid[i] )

import matplotlib.pyplot as plt
import matplotlib
plt.rcParams.update({
    "text.usetex": True,
    "font.size": 14}
)
plt.scatter(x, y, label="data")
plt.plot(xgrid, ygrid, label="interpolated")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("IMG_chapra_example_18_10.pdf")

