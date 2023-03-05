from scipy.optimize import linprog

c = [-150, -175]
A = [
    [7, 11],
    [10, 8]
]
b = [77, 80]

x0_bounds = (0, 9)
x1_bounds = (0, 6)

res = linprog(c, A_ub=A, b_ub=b, bounds=[x0_bounds, x1_bounds])

# We search for max, the result is negative the minimum we found from linprog
print("Max value: ", -res.fun)
print("Decision variables: ", res.x)
print("Message: ")
print(res.message)
