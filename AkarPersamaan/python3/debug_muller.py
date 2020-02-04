import numpy as np
import matplotlib.pyplot as plt

# Chapra 7th ed. halaman 200

def quadratic_poly(a, b, c, x2, x):
    return a*(x - x2)**2 + b*(x - x2) + c

def f(x):
    return x**3 - 13*x - 12

# Tebakan akar
x0 = 4.5
x1 = 5.5
x2 = 5.0

h0 = x1 - x0
h1 = x2 - x1

d0 = ( f(x1) - f(x0) ) / h0
d1 = ( f(x2) - f(x1) ) / h1

# Koefisien parabola
a = (d1 - d0) / (h1 + h0)
b = a*h1 + d1
c = f(x2)

print("a = ", a)
print("b = ", b)
print("c = ", c)

xmin = 0.0
xmax = 6.0
x_plot = np.linspace(xmin, xmax, 500)
fx = f(x_plot)

plt.clf()
plt.plot(x_plot, fx, label="f(x)")
# Tandai tebakan awal akar
plt.plot([x0], [f(x0)], marker="o", color="k")
plt.plot([x1], [f(x1)], marker="o", color="k")
plt.plot([x2], [f(x2)], marker="o", color="k")
#
plt.plot([x0], [0], marker="o", color="k")
plt.plot([x1], [0], marker="o", color="k")
plt.plot([x2], [0], marker="o", color="k")
# parabola
print("c = ", c)
print("f(x2) = ", f(x2))
print(quadratic_poly(a, b, c, x2, x0), " compare ", f(x0))
print(quadratic_poly(a, b, c, x2, x2), " compare ", f(x2))
f2x = quadratic_poly(a, b, c, x2, x_plot)
plt.plot(x_plot, f2x)
#
plt.ylim(-39, 120)
#
plt.grid()
plt.savefig("IMG_muller.png", dpi=150)
