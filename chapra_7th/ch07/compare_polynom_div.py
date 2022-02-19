import numpy as np
from numpy.polynomial.polynomial import polydiv

def my_polydiv(a, d):
    
    n = len(a) - 1
    m = len(d) - 1

    q = np.empty(n+1); q[:] = np.nan
    r = np.empty(n+1); r[:] = np.nan

    #print("initial: q = ", q)
    #print("initial: r = ", r)

    for j in range(n+1):
        r[j] = a[j]
        q[j] = 0.0
    
    for k in range(n-m,-1,-1): # step = -1, stop = 0 - 1
        q[k+1] = r[m+k]/d[m]
        for j in range(m+k-1,k-1,-1):
            r[j] = r[j] - q[k+1]*d[j-k]

    print("Before zeroing r: ", r)
    for j in range(m,n+1):
        r[j] = 0.0

    return q, r

#a = np.array([-41.0, 7.0, 1.0])
#d = np.array([-4.0, 3.0])

a = np.array([-41.0, 7.0, 1.0, 3.0])
d = np.array([-4.0, 3.0, 1.0])

print("a = ", a)
print("d = ", d)

# np.polydiv require reverse order
#q, r = np.polydiv(a[::-1], d[::-1])
#print("np.polydiv result: ")
#print("q = ", q[::-1])  # quotient
#print("r = ", r[::-1])  # remainder

# New API
q, r = polydiv(a, d)
print("polydiv result: ")
print("q = ", q)  # quotient
print("r = ", r)  # remainder

q, r = my_polydiv(a, d)
print("my_polydiv result: ")
print("a = ", a)
print("d = ", d)
print("q = ", q)  # quotient
print("r = ", r)  # remainder