from scipy.optimize import linprog


# minimize: c @ x
# such that:
# A_ub @ x <= b_ub
# A_eq @ x == b_eq
# lb <= x <= ub


# min -1*x0 + 4*x1
# such that
#                     -3*x0 +   x1 <= 6
# -x0 - 2*x1 >= -4 or    x0 + 2*x1 <= 4
# x1 >= -3
#
# bound of x0 is None
# Matrix A_eq is None
# Vector b_eq is None


c = [-1, 4]
A = [
    [-3, 1],
    [1, 2]
]
b = [6, 4]

x0_bounds = (None, None)
x1_bounds = (-3, None)

res = linprog(c, A_ub=A, b_ub=b, bounds=[x0_bounds, x1_bounds])
print("Min value: ", res.fun)
print("Decision variables: ", res.x)
print("Message: ")
print(res.message)
