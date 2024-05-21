# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# <h1 style="text-align: center">
# Persamaan Nilai Eigen
# </h1>

# %%
import numpy as np

# %%
A = np.array([
    [3.82, -1.75, 0],
    [-1.75, 3.82, -1.75],
    [0.0, -1.75, 3.82]
])

# %%
A

# %%
λ, X = np.linalg.eig(A) 

# %%
λ

# %%
X

# %%
X[:,0]

# %%
λ[0]

# %%
A @ X[:,0] 

# %%
λ[0] * X[:,0] 

# %%
A @ X[:,1] 

# %%
λ[1] * X[:,1] 

# %%
A @ X[:,2] 

# %%
λ[2] * X[:,2] 

# %%
np.dot( X[:,0], X[:,0] ) 

# %%
np.dot( X[:,0], X[:,1] ) 

# %%
np.dot( X[:,0], X[:,2] ) 

# %%
np.dot( X[:,1], X[:,1] ) 

# %%
np.dot( X[:,2], X[:,2] ) 

# %%
X.T @ X 


# %% [markdown]
# # Power method

# %%
def power_iteration(A, v0, NiterMax=1000, TOL=1e-8):
    v = np.copy(v0)
    for iterNo in range(NiterMax):
        v1 = A @ v
        # normalisasi
        v1[:] = v1[:]/np.linalg.norm(v1)
        dv = np.abs(v - v1)
        norm_dv = np.linalg.norm(dv)
        print("v1 = ", v1, "norm_dv = ", norm_dv)
        if norm_dv < TOL:
            print("Converged!")
            break
        v[:] = np.copy(v1)
    return v


# %%
v0 = np.array([1, 1, 1.0])

# %%
v = power_iteration(A, v0)

# %%
λ[0]

# %%
v.T @ ( A @ v )

# %%
v.T @ v 

# %%
X[:,0]

# %%
v2 = A @ v
v2

# %%
L = np.linalg.norm(v2)
L

# %%
v2n = v2[:]/L

# %%
np.linalg.norm(v2n)


# %%
def inverse_power_iteration(A, v0, NiterMax=1000, TOL=1e-8):
    v = np.copy(v0)
    Ainv = np.linalg.inv(A)
    for iterNo in range(NiterMax):
        v1 = Ainv @ v
        # normalisasi
        v1[:] = v1[:]/np.linalg.norm(v1)
        dv = np.abs(v - v1)
        norm_dv = np.linalg.norm(dv)
        print("v1 = ", v1, "norm_dv = ", norm_dv)
        if norm_dv < TOL:
            print("Converged!")
            break
        v[:] = np.copy(v1)
    return v


# %%
v = inverse_power_iteration(A, v0)

# %%
λ

# %%
v.T @ ( A @ v )

# %%
A

# %%
np.linalg.inv(A)

# %%
np.identity(3)

# %%
A.shape


# %%
def inverse_power_iteration(A, v0, NiterMax=1000, shift=0.0, TOL=1e-8):
    v = np.copy(v0)
    N = A.shape[0]
    I = np.identity(N)
    Ainv = np.linalg.inv(A - shift*I)
    λ = v.T @ (A @ v) / np.dot(v, v)
    for iterNo in range(NiterMax):
        v1 = Ainv @ v
        # normalisasi
        v1[:] = v1[:]/np.linalg.norm(v1)
        dv = np.abs(v - v1)
        norm_dv = np.linalg.norm(dv)
        λ_new = v.T @ (A @ v)
        delta_λ = abs(λ - λ_new)
        print("iterNo ", iterNo, " v1 = ", v1, " delta_λ = ", delta_λ)
        if delta_λ < TOL:
            print("Converged!")
            break
        v[:] = np.copy(v1)
        λ = λ_new
    return λ, v


# %%
λ, v = inverse_power_iteration(A, v0, shift=8)

# %%
λ

# %%
