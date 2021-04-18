import matplotlib.pyplot as plt
import numpy as np

SIGMA = 0.01

a = 1.0
b = 10.0

def my_func(x):
    return np.exp(-SIGMA*x)

def my_func2(t):
    x = 1/t
    return (1/t**2)*my_func(x)

t1 = 1/b
t2 = 1/a

xgrid = np.linspace(a, b, 1000)
tgrid = np.linspace(t1, t2, 1000)

plt.clf()
plt.plot(xgrid, my_func(xgrid))
plt.grid(True)
plt.savefig("IMG_transformed_v1_x.pdf")

plt.clf()
plt.plot(tgrid, my_func2(tgrid))
plt.grid(True)
plt.savefig("IMG_transformed_v1_t.pdf")