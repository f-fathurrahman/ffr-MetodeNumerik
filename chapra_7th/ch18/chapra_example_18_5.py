import numpy as np
from newton_interp import *

x = np.array([1.0, 4.0, 6.0 , 5.0, 3.0, 1.5, 2.5, 3.5])
y = np.array([0.0, 1.386294, 1.791759, 1.609438, 1.0986123, 0.4054641, 0.9162907, 1.2527630])

idx_choose = [0, 5, 6, 4]
print("idx_choose = ", idx_choose)

xi = 2.0

print("using the following data:")
print("x = ", x[idx_choose])
print("y = ", y[idx_choose])

yint, ea = newton_interp(x[idx_choose], y[idx_choose], xi)
  
print("yint = ", yint)
print("ea   = ", ea)

print("Error = ")
print(np.abs(yint[-1] - np.log(xi)))
