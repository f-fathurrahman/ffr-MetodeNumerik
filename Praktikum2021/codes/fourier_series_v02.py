import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
})

def sawtooth_wave(L, N, x):
    ss = 0.0
    for n in range(1,N+1):
        ss = ss + (1.0/n)*np.sin(n*np.pi*x/L)
    return 0.5 - ss/np.pi

N = 5
L = 1.0
NptsPlot = 500
x = np.linspace(0.0, 2*L, NptsPlot)
y = np.zeros(NptsPlot)
for i in range(NptsPlot):
    y[i] = sawtooth_wave(L, N, x[i])

plt.clf()
plt.plot(x, y)
plt.grid()
plt.savefig("IMG_fourier_series_v02.pdf")

