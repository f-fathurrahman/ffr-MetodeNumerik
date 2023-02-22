import numpy as np
import matplotlib.pyplot as plt

np.random.seed(13306023) # optional, for reproducibility

def my_func(x, y):
    return y - x - 2*x**2 - 2*x*y - y**2

xmin, xmax = -2.0, 2.0
ymin, ymax = 1.0, 3.0

xgrid = np.linspace(xmin, xmax, 200)
ygrid = np.linspace(ymin, ymax, 200)
X, Y = np.meshgrid(xgrid, ygrid)
FXY = my_func(X, Y)
plt.contour(X, Y, FXY, levels=10)
plt.gca().axis("equal")
plt.savefig("IMG_example_14_1.png", dpi=150)


Nsample = 10000

# First sample
x = xmin + (xmax - xmin)*np.random.rand()
y = ymin + (ymax - ymin)*np.random.rand()

fopt = my_func(x, y)
xopt, yopt = x, y

plt.plot([xopt], [yopt], marker="o", color="black")
plt.savefig("IMG_example_14_1_sample_0.png", dpi=150)

for i in range(1, Nsample+1):
    x = xmin + (xmax - xmin)*np.random.rand()
    y = ymin + (ymax - ymin)*np.random.rand()
    f = my_func(x, y)
    if f > fopt: # we search for max, if you search for min, use <
        fopt = f
        xopt, yopt = x, y
        plt.plot([xopt], [yopt], marker="o", color="black")
        plt.savefig("IMG_example_14_1_sample_" + str(i) + ".png", dpi=150)
    if i % 1000 == 0:
        print("Current no. sample = ", i)
        print("Current max: fopt = ", fopt, " at xopt = ", xopt, " yopt = ", yopt)


