# ---
# jupyter:
#   jupytext:
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

# %% [markdown]
# <h1 style="text-align:center;">TF2202 Komputasi Rekayasa - Persamaan Diferensial Parsial</h1>
# <h2 style="text-align:center;">Persamaan Laplace dan Poisson</h2>
# <h3 style="text-align:center;">Fadjar Fathurrahman</h3>

# %% [markdown]
# # Setup

# %%
import numpy as np
import matplotlib.pyplot as plt

# %%
# %matplotlib inline

# %%
from IPython.display import set_matplotlib_formats
set_matplotlib_formats("svg")

# %%
import matplotlib
matplotlib.style.use("default")

# %%
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

# %% [markdown]
# # Persamaan Laplace dan Poisson

# %% [markdown]
# Biasanya hanya bergantung pada variable spasial.

# %% [markdown]
# Persamaan Poisson:
# $$
# \nabla^2 u(x,y) = f(x,y)
# $$

# %% [markdown]
# Persamaan Laplace:
# $$
# \nabla^2 u(x,y) = 0
# $$

# %% [markdown]
# Untuk dua dimensi (spasial):
#
# $$
# \frac{\partial^2}{\partial x^2} u(x,y) +
# \frac{\partial^2}{\partial y^2} u(x,y)
# = f(x,y)
# $$

# %% [markdown]
# Gunakan notasi:
#
# - $u(x,y) = u_{i,j}$
#
# - $u(x+\Delta x, y) = u_{i+1,j}$
#
# - $u(x-\Delta x, y) = u_{i-1,j}$
#
# - $u(x, y + \Delta y) = u_{i,j+1}$
#
# - $u(x, y - \Delta y) = u_{i,j-1}$
#

# %% [markdown]
# Aproksimasi centered difference untuk turunan kedua:
#
# $$
# \frac{ u_{i+1,j} - 2u_{i,j} + u_{i-1,j} }{\Delta x^2} +
# \frac{ u_{i,j+1} - 2u_{i,j} + u_{i,j-1} }{\Delta y^2} =
# f_{i,j}
# $$

# %% [markdown]
# Persamaan ini dapat dituliskan menjadi sistem persamaan linear:
#
# $$
# \mathbf{A} \mathbf{u} = \mathbf{f}
# $$

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%
