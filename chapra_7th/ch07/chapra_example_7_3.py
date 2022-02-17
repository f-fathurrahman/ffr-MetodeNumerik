from root_bairstow import *

a = np.array([1.25, -3.875, 2.125, 2.75, -3.5, 1.0])

re, im = root_bairstow(a)

print("re = ", re)
print("im = ", im)