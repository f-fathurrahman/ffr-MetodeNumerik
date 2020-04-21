from math import pi

E = 10e9 # Pa
I = 1.25e-5 # m**4
L = 3 # m

# n = 1,2,3,...
def column_eigenvalue(n):
    return n*pi/L

# n = 1,2,3,...
def buckling_load(n):
    return n**2 * pi**2 * E * I / L**2

for n in range(1,9):
    print("%d %11.4f %11.3f" % (n, column_eigenvalue(n), buckling_load(n)/1e3))

