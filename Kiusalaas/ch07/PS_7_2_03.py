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

# %% [markdown]
# # Problem Set 7.1
#
# ## Soal 3

# %% [markdown]
# Selesaikan persamaan diferensial:
# $$
# y' = \sin(y)
# $$
# dengan syarat awal $y(0) = 1$ dari $x = 0$ sampai $x = 0.5$.
# Bandingkan dengan solusi analitik (fungsi implisit):
# $$
# x(y) = \mathrm{ln}\left( \mathrm{csc}(y) - \mathrm{cot}(y) \right) + 0.604582
# $$

# %%
import math


# %%
def csc(x):
    return 1/math.sin(x)


# %%
def cot(x):
    return 1/math.tan(x)


# %% [markdown]
# $$
# x(y) = \mathrm{ln}\left( \mathrm{csc}(y) - \mathrm{cot}(y) \right) + 0.604582
# $$

# %%
y = 1
math.log(csc(y) - cot(y)) + 0.604582

# %%
