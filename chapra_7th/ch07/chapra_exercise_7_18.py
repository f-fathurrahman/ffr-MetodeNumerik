from scipy.optimize import fsolve

def f(x_):
    # a = x_[0], u = x_[1], v = x_[2]
    a = x_[0]
    u = x_[1]
    v = x_[2]
    f1 = u**2 - 2*v**2 - a**2
    f2 = u + v - 2
    f3 = a**2 - 2*a - u
    return [f1, f2, f3]

root = fsolve(f, [0.0, 0.0, 0.0])
print(root)

# Test the value of function at root
print("f at root = ", f(root), " (should be close to zeros)")