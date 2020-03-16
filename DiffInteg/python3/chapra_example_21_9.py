import numpy as np

def my_func(x,y):
    return 2*x*y + 2*x - x**2 - 2*y**2 + 72

def integ_trapz( fa, fb, a, b ):
    I = (b - a) * (fa + fb) / 2
    return I

# 2 segment equal width
# h = (b - a)/2
def integ_simp13( f0, f1, f2, h ):
    I = h/3 * ( f0 + 4*f1 + f2 )
    return I

# Number of segments
Nx = 2
Ny = 2

xmin = 0.0
xmax = 8.0

ymin = 0.0
ymax = 6.0

area = (xmax - xmin) * (ymax - ymin)

x = np.linspace(xmin, xmax, Nx+1)
y = np.linspace(ymin, ymax, Ny+1)
hx = (xmax - xmin)/Nx
hy = (ymax - ymin)/Ny

I_y = np.zeros(Ny+1)
for j in range(Ny+1):
    I_y[j] = 0.0
    for i in range(Nx):
        fa = my_func(x[i], y[j])
        fb = my_func(x[i+1], y[j])
        I_y[j] = I_y[j] + integ_trapz(fa, fb, x[i], x[i+1])

I = 0.0
for j in range(Ny):
    I = I + integ_trapz( I_y[j], I_y[j+1], y[j], y[j+1] )

print()
print("Using trapezoidal rule: I   = %.10f" % I)
print("Using trapezoidal rule: avg = %.10f" % (I/area))

# Using Simpson 1/3 rule
I_y = np.zeros(Ny+1)
for j in range(Ny+1):
    # special case of Nx=2
    f0 = my_func(x[0], y[j])
    f1 = my_func(x[1], y[j])
    f2 = my_func(x[2], y[j])
    I_y[j] = integ_simp13( f0, f1, f2, hx )

I = integ_simp13( I_y[0], I_y[1], I_y[2], hy)

print()
print("Using Simpson 1/3 rule: I   = %.10f" % I)
print("Using Simpson 1/3 rule: avg = %.10f" % (I/area))

# Using SymPy
import sympy
x, y = sympy.symbols("x y")
f = 2*x*y + 2*x - x**2 - 2*y**2 + 72
I_sympy = sympy.integrate( f, (x, 0, 8), (y, 0, 6) )

print()
print("Using SymPy:            I   = %.10f" % I_sympy)
print("Using SymPy:            avg = %.10f" % (I_sympy/area))
