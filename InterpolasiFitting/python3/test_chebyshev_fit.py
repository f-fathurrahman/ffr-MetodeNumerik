import numpy as np
import matplotlib.pyplot as plt

from chebyshev_poly import *
from gen_chebyshev_nodes import *

z = gen_chebyshev_nodes(15)
print("z = ", z)
a = -1.0
b =  1.0

#x = np.array([0.75, 2, 3, 4, 6, 8, 8.5])
#y = np.array([1.2, 1.95, 2, 2.4, 2.4, 2.7, 2.6])

x = a + (b - a)*(z + 1)/2
#y = np.sin(2*np.pi*x/(b - a))
y = 1.0/(1 + 35*x**2)

print("x = ", x)
print("y = ", y)

Ndata = len(x)

plt.clf()
plt.plot(x, y, linewidth=0, marker="o", label="data")

NptsPlot = 200
A = np.min(x)
B = np.max(x)
x_plt = np.linspace(A,B,NptsPlot)
y_plt = np.zeros(NptsPlot)

for degree in range(1,Ndata+1,3):
    # fitting
    coefs = chebyshev_fit(x, y, degree)
    # plot
    for i in range(NptsPlot):
        y_plt[i] = chebyshev_fit_eval(coefs, A, B, x_plt[i])
    plt.plot(x_plt, y_plt, label="chebyshev-fit-"+str(degree))
    # hitung error
    err = 0.0
    for i in range(Ndata):
        y_approx = chebyshev_fit_eval(coefs, A, B, x[i])
        err = err + (y[i] - y_approx)**2
    print("degree = %d, err squared = %f" % (degree, err)) 

plt.legend()
plt.savefig("TEMP_test_chebyshev_fit.pdf")


