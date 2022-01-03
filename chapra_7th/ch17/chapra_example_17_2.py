import numpy as np

def fit_linear(x, y):
    N = len(x)
    num = N*np.sum(x * y) - np.sum(x) * np.sum(y)
    denum = N*np.sum(x**2) - np.sum(x)**2
    a1 = num/denum
    a0 = np.mean(y) - a1*np.mean(x)
    return a0, a1

data = np.loadtxt("table_17_1.dat")
x = data[:,0]
y = data[:,1]

N = len(x)
a0, a1 = fit_linear(x, y)

μ = np.mean(y)
St = np.sum((y - μ)**2)
sy = np.sqrt(St/(N-1))
print("sy = ", sy)

# Equation 17.8
Sr = np.sum( (y - a0 - a1*x)**2 )
print("Sr = ", Sr)

# Equation 17.9
syx = np.sqrt(Sr/(N-2))
print("syx = ", syx)

# Equation 17.10
r2 = (St - Sr)/St
print("r2 = ", r2)