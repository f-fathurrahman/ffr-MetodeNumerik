from sympy import *
x, y = symbols("x y")

eq1 = Eq( (x - 4)**2 + (y - 4)**4, 5 )
eq2 = Eq( x**2 + y**2, 16 )
plt1 = plot_implicit(eq1, show=False)
plt1.extend(plot_implicit(eq2, show=False))
plt1.save("IMG_chapra_exercise_7_13_sympy.pdf")

