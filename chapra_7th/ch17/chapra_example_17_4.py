import numpy as np
import matplotlib.pyplot as plt

def fit_linear(x, y):
    N = len(x)
    assert N == len(y)
    num = N*np.sum(x * y) - np.sum(x) * np.sum(y)
    denum = N*np.sum(x**2) - np.sum(x)**2
    a1 = num/denum
    a0 = np.mean(y) - a1*np.mean(x)
    return a0, a1

def calc_r2(x, y, a0, a1):
    N = len(x)
    assert N == len(y)
    μ = np.mean(y)
    St = np.sum((y - μ)**2)
    sy = np.sqrt(St/(N-1))
    Sr = np.sum( (y - a0 - a1*x)**2 )
    syx = np.sqrt(Sr/(N-2)) # not needed for r2
    r2 = (St - Sr)/St
    return r2

data = np.loadtxt("table_17_3.dat")

x = data[:,0]
y = data[:,1]

plt.clf()
plt.plot(np.log(x), np.log(y), marker="o", label="log model")
plt.grid(True)
plt.savefig("IMG_chapra_example_17_4.pdf")

xt = np.log10(x)
yt = np.log10(y)
a0, a1 = fit_linear(xt, yt)
print("a0 = ", a0)
print("a1 = ", a1)
print("r2 = ", calc_r2(xt, yt, a0, a1))