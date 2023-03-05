from scipy.optimize import linprog

c = [-6.0, -8.0]
A = [
    [5, 2],
    [6, 6],
    [2, 4]
]
b = [40, 60, 32]

x0_bounds = (0, None)
x1_bounds = (0, None)

res = linprog(c, A_ub=A, b_ub=b, bounds=[x0_bounds, x1_bounds])

# We search for max, the result is negative the minimum we found from linprog
print("Max value: ", -res.fun)
print("Decision variables: ", res.x)
print("Message: ")
print(res.message)
