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
# <h1 style="text-align:center;">TF2103 Komputasi Rekayasa - Persamaan Diferensial Parsial</h1>
# <h2 style="text-align:center;">Persamaan Laplace dan Poisson</h2>
# <h3 style="text-align:center;">Fadjar Fathurrahman</h3>

# %% [markdown]
# # Setup

# %%
import numpy as np
import matplotlib.pyplot as plt

# %%
# %matplotlib ipympl

# %%
import matplotlib_inline
matplotlib_inline.backend_inline.set_matplotlib_formats("svg")

# %%
import matplotlib
matplotlib.style.use("dark_background")
#matplotlib.rcParams.update({
#    "axes.grid" : True,
#    "grid.color": "gray"
#})

# %% [markdown]
# # Persamaan Laplace dan Poisson

# %% [markdown]
# Biasanya hanya bergantung pada variable spasial.

# %% [markdown]
# Persamaan Helmholtz (tipe PDE eliptik):
# $$
# \nabla^2 u(x,y) + g(x,y)u(x,y) = f(x,y)
# $$
# pada domain
# $$
# D = \left\{ (x,y) | x_0 \le x \le x_f, y_0 \le y \le y_f \right\}
# $$
# dengan kondisi batas Dirichlet (nilai fungsi diberikan pada titik-titik batas):
# $$
# \begin{align}
# u(x_0, y) & = b_{x_0}(y) \\
# u(x_f, y) & = b_{x_f}(y) \\
# u(x, y_0) & = b_{y_0}(x) \\
# u(x, y_f) & = b_{y_f}(x)
# \end{align}
# $$

# %% [markdown]
# Jika $g(x,y)=0$ dan $f(x,y) \neq 0$, diperoleh Persamaan Poisson:
# $$
# \nabla^2 u(x,y) = f(x,y)
# $$

# %% [markdown]
# Jika $g(x,y)=0$ dan $f(x,y) = 0$, diperoleh Persamaan Laplace:
# $$
# \nabla^2 u(x,y) = 0
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
# $$
# \frac{ u_{i+1,j} - 2u_{i,j} + u_{i-1,j} }{\Delta x^2} +
# \frac{ u_{i,j+1} - 2u_{i,j} + u_{i,j-1} }{\Delta y^2} +
# g_{i,j} u_{i,j}
# = f_{i,j}
# $$

# %% [markdown]
# Persamaan ini dapat dituliskan menjadi sistem persamaan linear yang dapat diselesaikan dengan metode standard (eliminasi Gauss, iteratif, dan sebagainya).
# $$
# \mathbf{A} \mathbf{u} = \mathbf{f}
# $$

# %% [markdown]
# Namun, kita akan menggunakan metode yang lebih sederhana, dengan skema iterasi sebagai berikut:
# $$
# u_{i,j}=\left[r_{y}\left(u_{i+1,j}+u_{i-1,j}\right)+r_{x}\left(u_{i,j+1}+u_{i,j-1}\right)+r_{xy}\left(g_{i,j}u_{i,j}-f_{i,j}\right)\right]
# $$
# dengan
# $$
# \begin{align}
# r_x & = \frac{\Delta x^{2}}{2\left(\Delta x^{2}+\Delta y^{2}\right)} \\
# r_y & = \frac{\Delta y^{2}}{2\left(\Delta x^{2}+\Delta y^{2}\right)} \\
# r_{xy} &= \frac{\Delta x^{2}\Delta y^{2}}{2\left(\Delta x^{2}+\Delta y^{2}\right)}
# \end{align}
# $$

# %% [markdown] jp-MarkdownHeadingCollapsed=true
# ## Implementasi

# %%
import numpy as np

def solve_poisson2d_dirichlet( \
    f_func, g_func, \
    bx0, bxf, by0, byf, \
    D, Nx, Ny, \
    TOL=1e-8, NiterMax=500, verbose=False):

    # Get boundary information
    x0 = D[0]; xf = D[1]
    y0 = D[2]; yf = D[3]

    Δx = (xf - x0)/Nx
    Δy = (yf - y0)/Ny

    x = np.linspace(x0, xf, Nx+1)
    y = np.linspace(y0, yf, Ny+1)

    u = np.zeros( (Nx+1, Ny+1) )

    for j in range(Ny+1):
        u[0,j] = bx0( y[j] )
        u[Nx,j] = bxf( y[j] )

    for i in range(Nx+1):
        u[i,0] = by0( x[i] )
        u[i,Ny] = byf( x[i] )


    # Initial values for other nodes
    sum_of_bv = np.sum(u[0,:]) + np.sum(u[Nx,:]) + \
                np.sum(u[:,0]) + np.sum(u[:,Ny])
    u[1:Nx,1:Ny] = sum_of_bv / (2*(Nx+Ny-2))

    f = np.zeros( (Nx+1, Ny+1) )
    g = np.zeros( (Nx+1, Ny+1) )
    for i in range(Nx+1):
        for j in range(Ny+1):
            f[i,j] = f_func(x[i], y[j])
            g[i,j] = g_func(x[i], y[j])

    rxy = 0.5 * Δx**2 * Δy**2 / (Δx**2 + Δy**2)
    ry = 0.5 * Δy**2 / (Δx**2 + Δy**2)
    rx = 0.5 * Δx**2 / (Δx**2 + Δy**2)
    for iterJacobi in range(NiterMax):
        err = 0.0
        u_old = np.copy(u)
        # loop only for internal nodes
        for j in range(1,Ny):
            for i in range(1,Nx):
                u[i,j] = rx*(u[i,j+1] + u[i,j-1]) + ry*(u[i+1,j] + u[i-1,j]) + \
                         rxy*( g[i,j]*u[i,j] - f[i,j])
        # calculate error
        err = np.mean( (u - u_old)**2 )
        if verbose:
            print("iterJacobi = %8d, err = %18.10e" % (iterJacobi, err))
        if err < TOL:
            print(f"Convergence is achieved in {iterJacobi} iterations, RMSE = {err}")
            break
    return u, x, y


# %% [markdown]
# ## Contoh 1

# %%
def f_func(x,y):
    return 0.0

def g_func(x,y):
    return 0.0

def bx0(y):
    return 0.0

def bxf(y):
    return 0.0

def by0(x):
    return np.sin(2*x)

def byf(x):
    return 0.0

Nx = 50
Ny = 50
D = [0.0, np.pi, 0.0, np.pi]

u, x, y = solve_poisson2d_dirichlet( \
    f_func, g_func, \
    bx0, bxf, by0, byf, \
    D, Nx, Ny, \
    TOL=1e-8, NiterMax=500)


# %%
X, Y = np.meshgrid(x, y)

fig = plt.figure()
plt.clf()
ax = fig.subplots(subplot_kw={"projection": "3d"})
plt.ylabel("y label")
plt.xlabel("x label")

# Need to transpose u
surf = ax.plot_surface(X, Y, u.T, linewidth=1, cmap=matplotlib.cm.coolwarm, antialiased=True)
#ax.view_init(elev=10, azim=135, roll=0)

# %%

# %%

# %%

# %%

# %%

# %%
