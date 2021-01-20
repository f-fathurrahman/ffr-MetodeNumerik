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
        ss = ss + (1.0/n)*np.sin(2*np.pi*n*x/L)
    return 0.5 - ss/np.pi


def square_wave(L, N, x):
    ss = 0.0
    for n in range(1,N+1,2):
        ss = ss + (1.0/n)*np.sin(2*np.pi*n*x/L)
    return 4*ss/np.pi



L = 1.0
NptsPlot = 1000

def do_plot(N, func):
    x = np.linspace(0.0, L, NptsPlot)
    y = np.zeros(NptsPlot)
    for i in range(NptsPlot):
        y[i] = func(L, N, x[i])
    plt.plot(x, y, label="N="+str(N))

plt.clf()
for n in [1,2,5,10,100]:
    #do_plot(n, sawtooth_wave)
    do_plot(n, square_wave)

plt.legend()
plt.grid()
plt.savefig("IMG_fourier_series_v03.pdf")