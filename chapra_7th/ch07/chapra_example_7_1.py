import numpy as np

# divide polynomial with factor (x-t)
def polynom_deflate(a_, t):
    a = np.copy(a_)
    #
    n = len(a) - 1 # degree
    r = a[n]
    a[n] = 0.0
    #
    for i in range(n-1,-1,-1): # loop until n-1 to 0
        s = a[i]
        a[i] = r
        print("i = ", i, end=" ")
        print("a = ", a)
        r = s + r*t
    return a, r

a = np.array([-40.0, 2.0, 1.0])

q, r = polynom_deflate(a, 4.0)
print("q = ", q)
print("r = ", r)