import numpy as np
import matplotlib.pyplot as plt

def my_func(x,y):
    return x*y**2

def grad_my_func(x,y):
    dfdx = y**2
    dfdy = 2*x*y
    return dfdx, dfdy


x = np.linspace(0.0, 4.0, 100)
y = np.linspace(0.0, 4.0, 100)
X, Y = np.meshgrid(x, y)
FXY = my_func(X, Y)

# Calc the gradient at (2,2)
x0, y0 = 2.0, 2.0
gx, gy = grad_my_func(x0, y0)

fig, ax = plt.subplots()
my_contour = ax.contour(X, Y, FXY, levels=np.linspace(8.0, 40.0, 5), colors="black")
ax.quiver(x0, y0, gx, gy, color="blue")
ax.set_aspect("equal", "box")
ax.clabel(my_contour, inline=True, fontsize=10)

plt.savefig("IMG_chapra_example_14_2.pdf")
plt.show()