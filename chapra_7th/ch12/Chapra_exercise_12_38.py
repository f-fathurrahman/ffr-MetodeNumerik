from sympy import *
x = symbols("x")
T = symbols("T", cls=Function)
h = symbols("h")
Ta = symbols("Ta")
sol = dsolve(
    T(x).diff(x,2) + h*( Ta - T(x) ),
    T(x),
    ics={T(0): 40.0, T(10): 200}
)
pprint(simplify(sol.subs({h: 0.02, Ta: 20.0})))