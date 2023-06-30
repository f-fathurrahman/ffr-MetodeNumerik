from sympy import *
init_printing()

b_1, b_2, b_3 = symbols("b_1 b_2 b_3", real=True)
f_1, f_2, f_3 = symbols("f_1 f_2 f_3", real=True)
x_1, x_2, x_3 = symbols("x_1 x_2 x_3", real=True)

a_2, a_3 = symbols("a_2 a_3", real=True)
c_1, c_2 = symbols("c_1 c_2", real=True)

eqn1 = Equality(b_1*x_1 + c_1*x_2, f_1)
eqn2 = Equality(a_2*x_1 + b_2*x_2 + c_2*x_3, f_2)
eqn3 = Equality(a_3*x_2 + b_3*x_3, f_3)

"""
gamma_1, gamma_2, gamma_3 = symbols("gamma_1 gamma_2 gamma_3")
beta_1, beta_2, beta_3 = symbols("beta_1 beta_2 beta_3")

# Recurrence relation
dict_subs = {
    x_1: gamma_1 - c_1/beta_1 * x_2
    x_2: gamma_2 - c_2/beta_2 * x_3
    x_3: gamma_3
}
"""

A = Matrix([
    [b_1, c_1, 0],
    [a_2, b_2, c_2],
    [0,   a_3, b_3], 
])
b = Matrix([
    [f_1],
    [f_2],
    [f_3]
])


# Row reduction
print()
α = A[1,0]/A[0,0]
A[1,:] = A[1,:] - A[0,:]*α
pprint(A)
b[1] = b[1] - b[0]*α
pprint(b)

print()
α = A[2,1]/A[1,1]
A[2,:] = A[2,:] - A[1,:]*α
pprint(A)
b[2] = b[2] - b[1]*α
pprint(b)
