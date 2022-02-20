from sympy import *

x = symbols("x")

a = []
Ndeg = 3
for i in range(0,Ndeg+1):
    a.append(symbols("a_" + str(i)))

p1 = a[0]
for i in range(1,Ndeg+1):
    p1 += a[i]*x**i
p1 = Poly(p1, x) # convert to a Polynomial object

t = symbols("t")
# divisor
d = x - t

b = []
Ndegq = Ndeg - 1
for i in range(0,Ndegq+1):
    b.append(symbols("b_" + str(i)))

q1 = b[0]
for i in range(1,Ndegq+1):
    q1 += b[i]*x**i

r = symbols("r") # remainder

p2 = expand(q1*d) + r
p2 = Poly(p2, x)

coefs1 = p1.all_coeffs()
coefs2 = p2.all_coeffs()

equations = []
# in reverse order, because coefs are given in reverse order (coef with high power first)
for i in range(Ndeg,-1,-1):
    equations.append(Equality(coefs1[i], coefs2[i]))

solb = [None for i in range(0,Ndegq+1)]

#In [92]: solve(equations[3], b[2])
#Out[92]: [a_3]

#In [91]: solve(equations[2], b[1])
#Out[91]: [a_2 + b_2*t]

#In [90]: solve(equations[1], b[0])
#Out[90]: [a_1 + b_1*t]

# Remainder
#In [89]: solve(equations[0], r)
#Out[89]: [a_0 + b_0*t]
