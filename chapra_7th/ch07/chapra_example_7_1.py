import numpy as np

# Algorithm to divide a polynomial (defi ned by its coefficients a)
# by a lower-order polynomial d.

def poly_division(a, n, d, m, q, r):
    
    for j in range(n+1):
        r[j] = a[j]
        q[j] = 0.0
    
    for k in range(n-m,-1,-1): # step = -1, stop = 0 - 1
        #print("k = ", k)
        q[k+1] = r[m+k]/d[m]
        for j in range(m+k-1,k-1,-1):
            r[j] = r[j] - q[k+1]*a[j-k]

    for j in range(m,n+1):
        r[j] = 0.0
    #print("r = ", r)
    #for i in range(0,n-m+1):
    #    a[i] = q[i+1]
    return

a = np.array([-24.0, 2.0, 1.0]) # x^2 + 2*x - 24
n = len(a) - 1

d = np.array([-4.0, 1.0]) # x - 4
m = len(d) - 1

q = np.zeros(n+1)
r = np.zeros(n+1)

print("a = ", a)
print("d = ", d)
print("q = ", q)
print("r = ", r)

poly_division(a, n, d, m, q, r)

print("after: ")
print("a = ", a)
print("d = ", d)
print("q = ", q)
print("r = ", r)

