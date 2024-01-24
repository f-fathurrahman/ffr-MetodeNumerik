import numpy as np
import matplotlib.pyplot as plt
from newton_interp import *

from chebyshev_poly import *
from gen_chebyshev_nodes import *

def func_01(x):
    return 1/(1 + 25*x**2)

a = -1.0
b =  1.0

Nnodes = 21
z = gen_chebyshev_nodes(Nnodes)
print("z = ", z)
x1 = a + (b - a)*(z + 1)/2
#x1 = np.linspace(a, b, Nnodes+1)
print("x1 = ", x1)
y1 = func_01(x1)

# dense grid
NptsPlot = 100
x2 = np.linspace(a, b, NptsPlot)

coefs = create_newton_polynom(x1, y1)
y2 = np.zeros((NptsPlot,))
for i in range(NptsPlot):
    y2[i] = eval_newton_polynom(coefs, x1, x2[i])

diff = y2 - func_01(x2)

plt.clf()
plt.plot(x1, y1, marker="o", label="sampled")
plt.plot(x2, y2, label="interp")
plt.plot(x2, func_01(x2), label="true")
plt.grid()
plt.legend()
plt.savefig("TEMP_newton_interp_02.pdf")

plt.clf()
plt.plot(x2, diff, marker="o", label="y2")
plt.grid()
plt.savefig("TEMP_newton_interp_02_diff.pdf")