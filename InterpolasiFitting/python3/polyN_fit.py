"""
c = polyFit(xData,yData,m).
Returns coefficients of the polynomial
p(x) = c[0] + c[1]x + c[2]xˆ2 +...+ c[m]xˆm
that fits the specified data in the least
squares sense.
sigma = stdDev(c,xData,yData).
Computes the std. deviation between p(x)
and the data.
"""
import numpy as np
import math

# m: derajat polinomial
def polyN_fit(xData, yData, m):
    A = np.zeros((m+1,m+1))
    b = np.zeros(m+1)
    s = np.zeros(2*m+1)
    for i in range(len(xData)):
        temp = yData[i]
        for j in range(m+1):
            b[j] = b[j] + temp
            temp = temp*xData[i]
        temp = 1.0
        for j in range(2*m+1):
            s[j] = s[j] + temp
            temp = temp*xData[i]
    for i in range(m+1):
        for j in range(m+1):
            A[i,j] = s[i+j]

    return np.linalg.solve(A,b)

def calc_std_dev_polyN(coefs, xData, yData):
    polyN = np.poly1d(coefs[::-1])

    N = len(xData) - 1
    M = len(coefs) - 1
    sigma = 0.0
    for i in range(N+1):
        valp = np.polyval(polyN, xData[i])
        sigma = sigma + (yData[i] - p)**2
    sigma = math.sqrt(sigma/(N-M))
    return sigma