from sympy import *

point1 = (1,2)
point2 = (2,1)
point3 = (4,5)

x1 = point1[0]; y1 = point1[1]
x2 = point2[0]; y2 = point2[1]
x3 = point3[0]; y3 = point3[1]

x = symbols("x")
fx = (x - y1)*(x - y3)*x1/(x1 - y1)*(x1 - y3) + \
     (x - y1)*(x - y2)*x2/(x2 - y2)*(x2 - y3) + \
     (x - y2)*(x - y3)*x3/(x3 - y2)*(x3 - y3)