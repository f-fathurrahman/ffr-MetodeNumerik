from sympy import *

x, t = symbols("x t")

f1 = 1/(1 + x**4)

f2 = f1.subs( {x: (1/t)**Rational(1,3)} )
pprint(f2)

my_integ = Integral( f1, (x,1,oo) )
pprint(my_integ)

#pprint(my_integ.doit())

# x**3 = 1/t  -> x = 1/x
# x -> 1 , t -> 1
# x -> oo, t -> 0

from mpmath import quad
f = lambda x: 1/(1 + x**4)

integ_old = 0.0
for pp in range(1,8):
    integ_res = quad(f, [1, 10**pp])
    print( pp, integ_res, abs(integ_res-integ_old) )
    integ_old = integ_res


f = lambda t: 1/(1 + (1/t)**(4/3))

print("Using change of variable: ", quad(f, [0.0, 1.0]) )
