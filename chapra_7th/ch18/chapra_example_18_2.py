import numpy as np
from newton_interp import *

x = np.array([1.0, 4.0, 6.0])
y = np.array([0.0, 1.386294, 1.791759])

xi = 2.0
  
yint, ea = newton_interp(x, y, xi)
  
print("yint = ", yint)
print("ea   = ", ea)
