import sympy as sym

x, t, a, L = sym.symbols("x t a L")

#u = x*(L-x)*5*t
u = x*2 * x*(L-x)*5*t

def pde(u):
    return sym.diff(u, t) - a*sym.diff(u, x, x)

f = sym.simplify(pde(u))
sym.pprint(f)
