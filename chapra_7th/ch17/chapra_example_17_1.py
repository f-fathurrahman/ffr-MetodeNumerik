import numpy as np

data = np.loadtxt("table_17_1.dat")
x = data[:,0]
y = data[:,1]
N = data.shape[0]

num = N*np.sum(x * y) - np.sum(x) * np.sum(y)
denum = N*np.sum(x**2) - np.sum(x)**2
a1 = num/denum
print("a1 = ", a1)

a0 = np.mean(y) - a1*np.mean(x)
print("a0 = ", a0)