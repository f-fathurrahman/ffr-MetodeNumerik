import numpy as np

def model_1(t):
    g = 9.8; m = 68.1; c = 12.5
    return g*m/c * (1 - np.exp(-c/m * t))

def model_2(t):
    g = 9.8; m = 68.1; c = 12.5
    return g*m/c * (t / (3.75 + t))

def fit_linear(x, y):
    N = len(x)
    num = N*np.sum(x * y) - np.sum(x) * np.sum(y)
    denum = N*np.sum(x**2) - np.sum(x)**2
    a1 = num/denum
    a0 = np.mean(y) - a1*np.mean(x)
    return a0, a1


dat = np.loadtxt("table_17_2.dat")
t_m = dat[:,0]
v_m = dat[:,1]

v_1 = model_1(t_m)
v_2 = model_2(t_m)

import matplotlib.pyplot as plt

plt.clf()
plt.plot(t_m, v_m, marker="o", label="measurement")
plt.plot(t_m, v_1, marker="^", label="model-1")
plt.plot(t_m, v_2, marker="s", label="model-2")
plt.xlabel("t")
plt.ylabel("v")
plt.legend()
plt.grid(True)
plt.savefig("IMG_chapra_example_17_3_t_v.pdf")

plt.clf()
plt.scatter(v_m, v_1, marker="o", label="model-1")
plt.scatter(v_m, v_2, marker="o", label="model-2")
plt.xlabel("v measurement")
plt.ylabel("v model")
plt.legend()
plt.grid(True)
ax = plt.gca()
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
plt.savefig("IMG_chapra_example_17_3.pdf")


model_1_a0, model_1_a1 = fit_linear(v_m, v_1)
model_2_a0, model_2_a1 = fit_linear(v_m, v_2)

print("Model 1: slope=%f intercept=%f" % (model_1_a1, model_1_a0))
print("Model 2: slope=%f intercept=%f" % (model_2_a1, model_2_a0))