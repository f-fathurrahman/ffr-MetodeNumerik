import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
})

L = 1.0
x = np.linspace(0,L,200)

suku1 = np.sin(2*np.pi*x/L)/1
suku2 = np.sin(2*3*np.pi*x/L)/3
suku3 = np.sin(2*5*np.pi*x/L)/5
suku4 = np.sin(2*7*np.pi*x/L)/7

plt.clf()
plt.plot(x, 4/np.pi*(suku1 + suku2 + suku3 + suku4))
plt.title(r"$f(x) = \frac{4}{\pi}\sum_{n=1,3,5,...}\frac{1}{n}\sin\left(\frac{n\pi x}{L}\right)$")
plt.grid(True)
plt.savefig("IMG_fourier_sqrwave.pdf")
