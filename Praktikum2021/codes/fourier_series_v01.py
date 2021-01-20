import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
})


# Eq.2 in Wikipedia
def fourier_series(a0, a, b, P, x):
    sN = a0/2
    N = len(a)
    # shoule check that len(a) == len(b)
    for n_ in range(N):
        n = n_ + 1
        t1 = a[n_]*np.cos(2*np.pi*n*x/P)
        t2 = b[n_]*np.sin(2*np.pi*n*x/P)
        sN = sN + t1 + t2
    return sN

A = 1.0
P = 2*np.pi # period

N = 100
a = np.zeros(N)
b = np.zeros(N)


# Fill up n
#a0 = 4*A/np.pi
#for n_ in range(N):
#    n = n_ + 1
#    is_even = n%2 == 0
#    is_odd = not is_even
#    if is_even:
#        a[n_] = -4*A/np.pi * ( 1.0/(n**2 - 1) )

#a0 = 2*A/np.pi
#for n_ in range(N):
#    n = n_ + 1
#    is_even = n%2 == 0
#    is_odd = not is_even
#    if is_even:
#        #a[n_] = -2*A/np.pi * ( 1.0/(1 - n**2) )
#        a[n_] = 2*A/np.pi * ( 1.0/(1 - n**2) )
#b[0] = A/2


#D = 0.3
#a0 = 2*A*D
#for n_ in range(N):
#    n = n_ + 1
#    is_even = n%2 == 0
#    a[n_] = A/(n*np.pi)*np.sin(2*np.pi*n*D)
#    b[n_] = 2*A/(n*np.pi)*(np.sin(np.pi*n*D))**2
#


#a0 = A
#for n_ in range(N):
#    n = n_ + 1
#    b[n_] = -A/(n*np.pi)

#a0 = A
#for n_ in range(N):
#    n = n_ + 1
#    b[n_] = A/(n*np.pi)


a0 = 2*A/3
for n_ in range(N):
    n = n_ + 1
    a[n_] = 4*A/(np.pi**2 * n**2)



NptsPlot = 400
x = np.linspace(0.0, 2*P, NptsPlot)  # plot two times the period
y = np.zeros(NptsPlot)
for i in range(NptsPlot):
    y[i] = fourier_series(a0, a, b, P, x[i])

plt.clf()
plt.plot(x, y)
plt.grid()
plt.savefig("IMG_fourier_series_v01.pdf")

print("Pass here")
