import numpy as np
import matplotlib.pyplot as plt

from numpy.random import RandomState
#rs = RandomState(1234)  # set random state for reproducibility
# OR
rs = RandomState()

# Global variables ...
SLOPE_EXACT = 2.4
INTERCEPT_EXACT = 0.6

def my_linear_func(x):
    return SLOPE_EXACT*x + INTERCEPT_EXACT

xmin = -5.0
xmax =  5.0
Nsamples = 20

x_sample = np.linspace(xmin, xmax, Nsamples)

# OR: (randomly spaced data)
#
#x_rand = (xmax - xmin)*rs.random_sample(Nsamples) + xmin
#x_sample = np.sort(x_rand)

y_exact = my_linear_func(x_sample)

NOISE_AMPLITUDE = 1.0
y_noisy = y_exact + rs.randn(Nsamples)*NOISE_AMPLITUDE

#plt.clf()
#plt.plot(x_sample, y_exact, label="exact", marker="o")
#plt.plot(x_sample, y_noisy, label="noisy", marker="o", linewidth=0)
#plt.legend()
#plt.savefig("IMG_DATA_v2.png", dpi=150)

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

scipy_res = scipy.stats.linregress(x_sample, y_noisy)

print("SLOPE_EXACT     = ", SLOPE_EXACT)
print("INTERCEPT_EXACT = ", INTERCEPT_EXACT)

print("\nSciPy result:")
print("slope     = ", scipy_res.slope)
print("intercept = ", scipy_res.intercept)
print("R-squared = ", scipy_res.rvalue**2)

 
x_plt = np.linspace(xmin, xmax, 200) # fine grid for plotting

plt.clf()
plt.plot(x_sample, y_noisy, label="data", marker="o", linewidth=0)
plt.plot(x_plt, scipy_res.slope*x_plt + scipy_res.intercept, label="fitted")
plt.legend()
plt.savefig("IMG_latihan_linreg_01.png", dpi=150)
