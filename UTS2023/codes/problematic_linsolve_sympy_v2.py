from sympy import *

#init_printing()

A00, A01, A10, A11 = symbols("A00 A01 A10 A11")
b0, b1 = symbols("b0 b1")

A = Matrix([
    [A00, A01],
    [A10, A11]
])

b = Matrix([[b0], [b1]])

x = A.solve(b)
print("x = ");
print(x)

print("Residual = ")
print(simplify(A * x - b))
