from gen_chebyshev_nodes import *

z = gen_chebyshev_nodes(5)
print("z = ", z)
a = -3.0
b =  2.0
x = a + (b - a)*(z + 1)/2
#x = (a + b - 2*z)/(a - b)
print("x = ", x)