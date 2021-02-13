import sympy
a, u, v = sympy.symbols("a u v")

eq1 = u**2 - 2*v**2 - a**2
eq2 = u + v - 2.0 # using real value to force numerical evaluation
eq3 = a**2 - 2*a - u

sols = sympy.nonlinsolve([eq1, eq2, eq3], [a, u, v])
# sols is finite set, we will iterate over its elements
for sol in sols:
    print("sol = ", sol)
