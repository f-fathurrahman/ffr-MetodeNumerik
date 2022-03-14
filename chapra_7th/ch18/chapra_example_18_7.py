import numpy as np
from lagrange_interp import *

t = np.array([ 1.0, 3.0, 5.0, 7.0, 13.0 ])
v = np.array([ 800.0, 2310.0, 3090.0, 3940.0, 4755.0 ])
tt = 10.0
vv = lagrange_interp(t, v, tt )

print("vv = ", vv)

