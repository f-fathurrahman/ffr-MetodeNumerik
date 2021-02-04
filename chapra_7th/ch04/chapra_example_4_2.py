from sympy import *

init_printing(use_unicode=True)
# if you are using Jupyter Lab or Notebook, use the following line instead:
#init_printing(use_latex=True)

x = symbols("x")

f = cos(x)

xi = pi/4
xip1 = pi/3
h = xip1 - xi

# zeroth order
f_approx = diff(f, x, 0).subs({x: xi}) # or simply call cos(xi)

for n in range(1,7): # from 1 to 6
    new_term = diff(f, x, n) * h**n / factorial(n)
    f_approx = f_approx + new_term
    pprint(f_approx)
    print(N(f_approx.subs({x: xi}))) # use N to force numerical expression

f_true = N(f.subs({x: xip1}))
print("f_true = ", f_true)
print(type(f_true))