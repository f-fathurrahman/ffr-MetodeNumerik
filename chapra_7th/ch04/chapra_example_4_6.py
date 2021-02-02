from sympy import *

F, L, E, I = symbols("F L E I")
y = F*L**4/(8*E*I)

F_num = 750; ΔF = 30 # N/m
L_num = 9; ΔL = 0.03 # m
E_num = 7.5e9; ΔE = 5e7 # N/m^2
I_num = 0.0005; ΔI = 0.000005 # m^4

dict_subs = {F: F_num, L: L_num, E: E_num, I: I_num}
y_num = y.subs(dict_subs)
print("y = ", y_num)

Δy = abs(diff(y,F))*ΔF + abs(diff(y,L))*ΔL + abs(diff(y,E))*ΔE + abs(diff(y,I))*ΔI

pprint(Δy)
print("Δy = ", Δy.subs(dict_subs))

