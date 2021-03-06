from sympy import *

init_printing(use_unicode=True)
# if you are using Jupyter Lab or Notebook, use the following line instead:
#init_printing(use_latex=True)

x, h = symbols("x h")

f = log(x)

xi = 1
xip1 = 2.5 #Rational(5,2)
h_num = xip1 - xi

# zeroth order
f_approx = diff(f, x, 0).subs({x: xi})
f_true = N(f.subs({x: xip1}))

dict_subs = {h: h_num, x: xi}

for n in range(1,10):
    new_term = diff(f, x, n) * h**n / factorial(n)
    f_approx = f_approx + new_term
    #pprint(f_approx)
    f_approx_num = N(f_approx.subs(dict_subs))
    #print(f_approx_num)
    pprint(new_term)
    new_term_num = N(new_term.subs(dict_subs))
    err = abs(f_true - f_approx_num)
    print("%18.10f %18.10f %18.10f" % (f_approx_num, new_term_num, err))
    #print("err = ", err)

print("h = ", h_num)
print("f_true = ", f_true)
print(type(f_true))