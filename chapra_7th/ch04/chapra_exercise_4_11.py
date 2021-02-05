from sympy import *

g = 9.81
m = 50
t = 6
c = symbols("c")
v = g*m/c*( 1 - exp(-(c/m)*t) )

c_num = 12.5
Δc = 1.5
dict_subs = {c: c_num}
v_num = v.subs(dict_subs)
print("v = ", v_num)

Δv = abs(diff(v,c))*Δc

pprint(Δv)
print("Δv = ", Δv.subs(dict_subs))

dict_subs = {
    c: c_num + Δc
}
vmin = v.subs(dict_subs)
print("vmin = ", vmin)
print("v - Δv = ", v_num - Δv.subs({c: c_num}))

dict_subs = {
    c: c_num - Δc
}
vmax = v.subs(dict_subs)
print("vmax = ", vmax)
print("v + Δv = ", v_num + Δv.subs({c: c_num}))
