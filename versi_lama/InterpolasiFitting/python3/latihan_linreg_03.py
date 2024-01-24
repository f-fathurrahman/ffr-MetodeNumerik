import numpy as np
import matplotlib.pyplot as plt

from numpy.random import RandomState
#rs = RandomState(1234)  # set random state for reproducibility
# OR
rs = RandomState()

# Global variables ...
α_EXACT = 2.4
β_EXACT = 1.6

def my_power_func(x):
    return α_EXACT* x**β_EXACT

xmin = 1e-2 # don't start from 0
xmax = 5.0
Nsamples = 20

#x_sample = np.linspace(xmin, xmax, Nsamples)
x_sample_transf = np.linspace(np.log(xmin), np.log(xmax), Nsamples)
x_sample = np.exp(x_sample_transf)

# OR: (randomly spaced data)
#
#x_rand = (xmax - xmin)*rs.random_sample(Nsamples) + xmin
#x_sample = np.sort(x_rand)

y_exact = my_power_func(x_sample)

NOISE_AMPLITUDE = 1e-2
y_noisy = y_exact + rs.randn(Nsamples)*NOISE_AMPLITUDE
# Additional check such that no negative value are encountered
idx = (y_noisy < 0.0)
y_noisy[idx] = np.abs(y_noisy[idx])

plt.clf()
plt.plot(x_sample, y_exact, label="exact", marker="o")
plt.plot(x_sample, y_noisy, label="noisy", marker="o", linewidth=0)
plt.legend()
plt.savefig("IMG_DATA_v3.png", dpi=150)

# Do linear regression here ...

# FILL IN YOUR CODE HERE ....
#def my_linear_regression(x, y):
#
#    # BLAH BLAH BLAH
#
#    # Return important results here ...
#    return slope, intercept, rsquared



# Using scipy.stats (for comparison)
# You should use your own function
import scipy.stats

# Transform the data to linear
x_transf = np.log(x_sample)
y_transf = np.log(y_noisy)

# linear regression
scipy_res = scipy.stats.linregress(x_transf, y_transf)

print("α_EXACT = ", α_EXACT)
print("β_EXACT = ", β_EXACT)


α = np.exp(scipy_res.intercept)
print("scipy_res.intercept = ", scipy_res.intercept)
print("log(α) = ", np.log(α))

β = scipy_res.slope
print("\nSciPy result:")
print("α = ", α)
print("β = ", β)
print("R-squared = ", scipy_res.rvalue**2)

 
x_plt = np.linspace(xmin, xmax, 200) # fine grid for plotting

plt.clf()
plt.plot(x_sample, y_noisy, label="data", marker="o", linewidth=0)
plt.plot(x_plt, α * x_plt**β, label="fitted")
plt.legend()
plt.savefig("IMG_latihan_linreg_03.png", dpi=150)

plt.clf()
plt.plot(x_transf, y_transf, marker="o", label="transformed")
x_plt_transf = np.log(x_plt)
plt.plot(x_plt_transf, scipy_res.intercept + scipy_res.slope*x_plt_transf, label="transformed fit")
plt.savefig("IMG_DATA_v3_transf.png", dpi=150)
