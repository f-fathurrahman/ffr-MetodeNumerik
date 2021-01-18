from scipy.optimize import fsolve

def f(x_):
    # x = x_[0], y = x_[1]
    x = x_[0]
    y = x_[1]
    f1 = (x - 4)**2 + (y - 4)**4 - 5
    f2 = x**2 + y**2 - 16
    return [f1, f2]

# First root
root = fsolve(f, [2.0, 4.0])
print(root)
print("f at root = ", f(root), " (should be close to zeros)")

# Second root
root = fsolve(f, [3.0, 2.0])
print(root)
print("f at root = ", f(root), " (should be close to zeros)")