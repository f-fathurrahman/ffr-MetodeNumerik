from scipy.optimize import fsolve

def f(x_):
    # x = x_[0], y = x_[1]
    x = x_[0]
    y = x_[1]
    f1 = -x**2 + x + 0.75 - y
    f2 = y + 5*x*y - x**2
    return [f1, f2]

root = fsolve(f, [1.2, 1.2])
print(root)

# Test the value of function at root
print("f at root = ", f(root), " (should be close to zeros)")