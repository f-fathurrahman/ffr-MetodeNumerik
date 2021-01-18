from scipy.optimize import fsolve

def f(x_):
    # x = x_[0] and y = x_[1] 
    x = x_[0]
    y = x_[1]
    u = x**2 + x*y - 10
    v = y + 3*x*y**2 - 57
    return [u, v]

#root = fsolve(f, [1.0, 1.0])
root = fsolve(f, [1.0, 3.5]) # using initial guess as in the book
print(root)
