import numpy as np

np.random.seed(13306023) # optional, for reproducibility

def my_func(x, y):
    return y - x - 2*x**2 - 2*x*y - y**2


xmin, xmax = -2.0, 2.0
ymin, ymax = 1.0, 3.0

Nsample = 10000

# First sample
x = xmin + (xmax - xmin)*np.random.rand()
y = ymin + (ymax - ymin)*np.random.rand()

fopt = my_func(x, y)
xopt, yopt = x, y

for i in range(1, Nsample+1):
    x = xmin + (xmax - xmin)*np.random.rand()
    y = ymin + (ymax - ymin)*np.random.rand()
    f = my_func(x, y)
    if f > fopt: # we search for max, if you search for min, use <
        fopt = f
        xopt, yopt = x, y
    if i % 1000 == 0:
        print("Current no. sample = ", i)
        print("Current max: fopt = ", fopt, " at xopt = ", xopt, " yopt = ", yopt)


