from sympy import Matrix, pprint, init_printing

init_printing()

A = Matrix([
    [888445, 887112],
    [887112, 885781]
])

b = Matrix([[10], [0]])

x = A.solve(b)
print("x = ");
pprint(x)

print("Residual = ")
pprint(A * x - b)
