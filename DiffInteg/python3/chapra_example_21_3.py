from math import exp

def my_func(t):
    g = 9.8
    m = 68.1
    c = 12.5
    return g*m/c * (1 - exp(-(c/m)*t) )

# N is number of segments
# Npoints is N + 1
def integ_trapz_multiple( f, a, b, N ):
    x0 = a
    xN = b
    h = (b - a)/N
    ss = 0.0
    for i in range(1,N):
        xi = x0 + i*h
        ss = ss + f(xi)
    I = h/2 * ( f(x0) + 2*ss + f(xN) )
    return I

t = 10.0
a = 0.0
b = t

d_exact = 289.43515

print("------------------------------------------------------")
print("   N          h          d          E_t         ε_t")
print("------------------------------------------------------")
for Nsegments in [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 5000, 10000]:
    h = (b - a)/Nsegments
    d = integ_trapz_multiple( my_func, a, b, Nsegments )
    E_t = d_exact - d
    ε_t = E_t/d_exact * 100
    print("%5d  %10.4f  %10.4f   %10.4e %10.2e%%" % (Nsegments, h, d, E_t, ε_t))


# Calculate "Exact" result using SymPy
import sympy
g = 9.8
m = 68.1
c = 12.5
t = sympy.symbols("t")
f = g*m/c * (1 - sympy.exp(-(c/m)*t) )
d_sympy = sympy.integrate( f, (t,0,10))

print()
print("SymPy result:")
print("d_sympy = ", d_sympy)
