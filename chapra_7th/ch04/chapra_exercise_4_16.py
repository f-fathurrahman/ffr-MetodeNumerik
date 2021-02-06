from sympy import *

B, H, S, n = symbols("B H S n")
# B = width, H = depth, S = slope, n = roughness

Q = 1/n * (B*H)**(5/3) / (B + 2*H)**(2/3) * sqrt(S)

B_num = 20
H_num = 0.3
n_num = 0.03; Δn = n_num*0.1
S_num = 0.0003; ΔS = S_num*0.1
dict_subs = {B: B_num, H: H_num, n: n_num, S: S_num}

Q_num = Q.subs(dict_subs)
print(Q_num)

ΔQ_n = abs(diff(Q,n).subs(dict_subs))*Δn
print("ΔQ_n = ", ΔQ_n)

ΔQ_s = abs(diff(Q,S).subs(dict_subs))*ΔS
print("ΔQ_s = ", ΔQ_s)