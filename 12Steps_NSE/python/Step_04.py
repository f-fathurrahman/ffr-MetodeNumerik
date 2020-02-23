import sympy

x, nu, t = sympy.symbols("x nu t")
phi = sympy.exp( -(x - 4*t)**2 / (4*nu*(t + 1)) ) + \
      sympy.exp( -(x - 4*t - 2*sympy.pi)**2 / (4 * nu * (t + 1)) )
sympy.pprint(phi)

phiprime = phi.diff(x)
sympy.pprint(phiprime)

u = -2*nu*(phiprime/phi) + 4
sympy.pprint(u)

from sympy.utilities.lambdify import lambdify
ufunc = lambdify((t,x,nu), u)
print(ufunc(1,4,3))

import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.style.use("dark_background")

Nx = 101
Nt = 100
dx = 2*np.pi/(Nx-1)
nu = 0.07
dt = dx*nu

x = np.linspace(0.0, 2*np.pi, Nx)
un = np.empty(Nx)
t = 0.0
u = np.asarray([ufunc(t, x0, nu) for x0 in x])

plt.clf()
plt.plot(x, u)
plt.xlim(0.0, 2*np.pi)
plt.ylim(0.0, 10.0)
plt.savefig("IMG_Step_04_0000.png", dpi=150)

for n in range(Nt):
    un = u.copy()
    for i in range(1, Nx-1):
        u[i] = un[i] - un[i] * dt / dx *(un[i] - un[i-1]) + nu * dt / dx**2 * (un[i+1] - 2 * un[i] + un[i-1])
        u[0] = un[0] - un[0] * dt / dx * (un[0] - un[-2]) + nu * dt / dx**2 * (un[1] - 2 * un[0] + un[-2])
        u[-1] = u[0]
    u_analytical = np.asarray([ufunc((n+1) * dt, xi, nu) for xi in x])
    plt.clf()
    plt.plot(x, u, label="numeric")
    plt.plot(x, u_analytical, label="analytic")
    plt.xlim(0.0, 2*np.pi)
    plt.ylim(0.0, 10.0)
    plt.legend()
    filename = "IMG_Step_04_{:04d}.png".format(n + 1)
    plt.savefig(filename, dpi=150)

