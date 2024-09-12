# %%
import numpy as np
import matplotlib.pyplot as plt

# %%
import matplotlib
matplotlib.style.use("dark_background")

# %%
import matplotlib_inline
matplotlib_inline.backend_inline.set_matplotlib_formats("svg")

# %% [markdown]
# Grid spasial:

# %%
x0 = 0.0 # ujung kiri
xL = 2.0 # ujung kanan
Nx = 41
Δx = (xL - x0)/(Nx - 1)
x = np.linspace(x0, xL, Nx)

# %% [markdown]
# Pastikan bahwa `Δx` merupakan jarak antara titik pada `x`.

# %% [markdown]
# Grid temporal (waktu)

# %%
Nt = 25
Δt = 0.025
c = 1.0  # kecepatan


# %%

# %%
def setup_initial_cond(xgrid):
    Nx = len(xgrid)
    u = np.ones(Nx)
    # u = 2 between 0.5 and 1 as per IC
    idx = (x > 0.5) & (x < 1)
    u[idx] = 2.0
    return u


# %%
u = setup_initial_cond(x)

# %%
plt.plot(x, u);

# %% [markdown]
# Baru kita implementasikan skema numerik:

# %%
u = setup_initial_cond(x)
un = np.zeros(Nx) # for previous time grid
# u is for the next
t = 0.0
for n in range(Nt):
    un = u.copy() # un is previous time grid
    for i in range(1,Nx):
        u[i] = un[i] - c*Δt/Δx * ( un[i] - un[i-1] )
    plt.clf()
    plt.plot(x, u)
    plt.ylim(0.9, 2.1)
    filename = "IMG_Step_01_{:04d}.png".format(n)
    t += Δt # increment time
    plt.title("t = {:.2f}".format(t))
    plt.savefig(filename, dpi=150)


# %%

# %%
