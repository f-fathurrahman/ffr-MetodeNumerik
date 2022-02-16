from numpy.polynomial import Polynomial

# coefficients are a0, a1, a2, ..., aN
# N is degree of the polinomial
p = Polynomial([1.25, -3.875, 2.125, 2.75, -3.5, 1.0])
print("p = ", p)
print(p.roots())