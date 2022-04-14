from sympy import *

x = symbols("x")
f = sin(x)
a = 0
b = 1
res = integrate(f, (x,0,1))

pprint(res)
print(N(res))
