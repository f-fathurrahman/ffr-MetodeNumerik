# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import numpy as np

# %%
import matplotlib.style
matplotlib.style.use("dark_background")

# %%
import matplotlib.pyplot as plt

# %% [markdown]
# Exact solution:
# $$
# u(x,t) = 5 t x(L-x)
# $$

# %%
α = 0.5 # koefisien difusi
L = 1.5 # panjang domain


# %%
def u_exact(x,t):
    return 5*t*x*(L - x)


# %%
def initial_cond(x):
    return u_exact(x, 0)


# %%
def boundary_cond_x0(t):
    return u_exact(0,t)


# %%
def boundary_cond_xL(t):
    return u_exact(L,t)


# %% [markdown]
# $$
# f(x,t) = 10\alpha t+5x(L-x)
# $$

# %%
def f_source(x,t):
    return 10*α*t + 5*x*(L - x)


# %%
Nx = 10
t_final = 0.1
Nt = 5
x = np.linspace(0, L, Nx + 1)
t = np.linspace(0, t_final, Nt + 1)

# %% [markdown]
# Fourier mesh number
# $$
# F=\frac{\alpha\Delta t}{\Delta x^{2}}
# $$

# %%
Δt = t[1] - t[0]
Δx = x[1] - x[0]

# %%
F = α * Δt / Δx**2
F

# %% [markdown]
# Metode forward Euler

# %% [markdown]
# $$
# u_{i}^{n+1}=u_{i}^{n}+F\left(u_{i+1}^{n}-2u_{i}^{n}+u_{i-1}^{n}\right)+\Delta t\,f_{i}^{n}
# $$

# %%
# Solusi
u = np.zeros( (Nx+1, Nt+1) )

# %%
# Set syarat awal
u[:,0] = initial_cond(x[:])

# %%
# Set syarat batas
u[0,0] = boundary_cond_x0(t[0])
u[Nx,0] = boundary_cond_xL(t[0])

# %%
n = 0
for i in range(1,Nx):
    u[i,n+1] = u[i,n] + F*( u[i+1,n] - 2*u[i,n] + u[i-1,n] ) + Δt * f_source(x[i], t[n])

# %%
plt.plot(x, u[:,1], label=f"t = {t[n+1]} numerik")
plt.plot(x, u_exact(x, t[n+1]), label=f"t = {t[n+1]} eksak")
plt.legend();

# %%
n = 1
for i in range(1,Nx):
    u[i,n+1] = u[i,n] + F*( u[i+1,n] - 2*u[i,n] + u[i-1,n] ) + Δt * f_source(x[i], t[n])

# %%
plt.plot(x, u[:,n+1], label=f"t = {t[n+1]} numerik")
plt.plot(x, u_exact(x, t[n+1]), label=f"t = {t[n+1]} eksak")
plt.legend();

# %%
n = 2
for i in range(1,Nx):
    u[i,n+1] = u[i,n] + F*( u[i+1,n] - 2*u[i,n] + u[i-1,n] ) + Δt * f_source(x[i], t[n])

# %%
plt.plot(x, u[:,n+1], label=f"t = {t[n+1]} numerik")
plt.plot(x, u_exact(x, t[n+1]), label=f"t = {t[n+1]} eksak")
plt.legend();

# %%
