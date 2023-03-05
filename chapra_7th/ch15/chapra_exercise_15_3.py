from scipy.optimize import linprog

c = [-1.75, -1.25]
A = [
    [1.2, 2.25],
    [1.0, 1.1],
    [2.5, 1]
]
b = [14, 8, 9]

x0_bounds = (0, None)
x1_bounds = (0, None)

res = linprog(c, A_ub=A, b_ub=b, bounds=[x0_bounds, x1_bounds])

# We search for max, the result is negative the minimum we found from linprog
print("Max value: ", -res.fun)
print("Decision variables: ", res.x)
print("Message: ")
print(res.message)
