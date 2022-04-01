import numpy as np
from newton_interp import *

x = np.array([1.0, 4.0, 6.0 , 5.0, 3.0, 1.5, 2.5, 3.5])
y = np.array([0.0, 1.386294, 1.791759, 1.609438, 1.0986123,
    0.4054641, 0.9162907, 1.2527630])

N = len(x)
xi = 2.0
yint, ea = newton_interp(x[::-1], y[::-1], xi)
  
true_val = np.log(xi)
print("True value = ", true_val)

print("-----------------------------------------------")
print("Order        yint          ea        true error")
print("-----------------------------------------------")
for i in range(0,N):
    if i != (N-1):
        print("%2d %18.10f %12.5e %12.5e" % (i, yint[i], ea[i], true_val-yint[i]))
    else:
        print("%2d %18.10f %12.5e %12.5e" % (i, yint[i], np.nan, true_val-yint[i]))
# We use if-else here because ea is only given up to N-1 for order N interpolation.