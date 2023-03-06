from scipy.optimize import minimize

def obj_func(x):
    x1 = x[0]
    x2 = x[1]
    # we use the same notation as in the book
    return 2 + x1 - x2 + 2*x1**2 + 2*x1*x2 + x2**2

x0 = [-0.5, 0.5] # initial search point

# The following are some methods that do not require gradient information
# result = minimize(obj_func, x0, method="Nelder-Mead")
# result = minimize(obj_func, x0, method="BFGS")
result = minimize(obj_func, x0, method="L-BFGS-B")

# NOTE
# We search for the minimum value.
# The book mistakenly asks us to find the maximum instead.

# Full output (might differ for each methods)
print("result = ")
print(result)

# The result
print("result.x (arg) = ", result.x)
print("result.fun (minimum value) = ", result.fun)
