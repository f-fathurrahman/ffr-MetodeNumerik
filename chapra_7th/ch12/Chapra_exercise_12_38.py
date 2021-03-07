from sympy import *
x = symbols("x")
T = symbols("T", cls=Function)
h = symbols("h") # 0.02
Ta = symbols("Ta") # 20
#h = 0.02
#Ta = 20.0
sol = dsolve(
    T(x).diff(x,2) + h*( Ta - T(x) ),
    T(x),
    #ics={T(0): 40.0, T(10): 200}
)
pprint(sol)