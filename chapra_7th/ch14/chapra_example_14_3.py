from sympy import *

x, y, h = symbols("x y h")

fxy = 2*x*y + 2*x - x**2 - 2*y**2

x0 = -1; y0 = 1
dict_subs = { x: x0, y: y0 }
dfdx = diff(fxy, x).subs(dict_subs)
dfdy = diff(fxy, y).subs(dict_subs)

dict_subs = {
    x: x0 + dfdx*h,
    y: y0 + dfdy*h
}
fh = fxy.subs(dict_subs)
pprint(simplify(fh))

