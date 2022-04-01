import numpy as np
from newton_interp import *

x = np.array([1.0, 4.0, 6.0 , 5.0, 3.0, 1.5, 2.5, 3.5])
y = np.array([0.0, 1.386294, 1.791759, 1.609438, 1.0986123,
    0.4054641, 0.9162907, 1.2527630])

N = len(x)
xi = 2.0
true_val = np.log(xi)
print("True value = ", true_val)

yint, ea1 = newton_interp(x, y, xi)
print()
print("Original ordering:")
print("-----------------------------------------------")
print("Order        yint          ea        true error")
print("-----------------------------------------------")
errors1 = true_val - yint # true error
for i in range(0,N):
    if i != (N-1):
        print("%2d %18.10f %12.5e %12.5e" % (i, yint[i], ea1[i], errors1[i]))
    else:
        print("%2d %18.10f %12.5e %12.5e" % (i, yint[i], np.nan, errors1[i]))
# We use if-else here because ea is only given up to N-1
# for order N interpolation.


yint, ea2 = newton_interp(x[::-1], y[::-1], xi)
print()
print("Reversed ordering:")
print("-----------------------------------------------")
print("Order        yint          ea        true error")
print("-----------------------------------------------")
errors2 = true_val - yint # true error
for i in range(0,N):
    if i != (N-1):
        print("%2d %18.10f %12.5e %12.5e" % (i, yint[i], ea2[i], errors2[i]))
    else:
        print("%2d %18.10f %12.5e %12.5e" % (i, yint[i], np.nan, errors2[i]))

import matplotlib.pyplot as plt
import matplotlib
plt.rcParams.update({
    "text.usetex": True,
    "font.size": 14}
)
plt.clf()
plt.plot(ea1, label="est. error (original)", marker="o")
plt.plot(ea2, label="est. error (reversed)", marker="o")
plt.plot(np.abs(errors1), label="abs. true error (original)", marker="o")
plt.grid(True)
plt.xlabel("Order of interpolation")
plt.ylabel("Error")
plt.legend()
plt.tight_layout()
plt.savefig("IMG_chapra_example_18_5.pdf")