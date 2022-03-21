import numpy as np
from newton_interp import *

x = np.array([ 1.0, 4.0, 6.0 , 5.0 ])
y = np.array([ 0.0, 1.386294, 1.791759, 1.609438 ])

xi = 2.0
  
yint, ea = newton_interp(x, y, xi)
  
print("yint  = %18.10f" % yint[-1])
print("ytrue = %18.10f" % np.log(xi))

err = np.abs(yint[-1] - np.log(xi))
print("err   = %18.10e" % err)
print("ea    = ", ea)