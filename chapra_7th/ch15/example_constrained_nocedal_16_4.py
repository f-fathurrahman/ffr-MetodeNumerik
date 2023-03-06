from scipy.optimize import minimize


fun = lambda x: (x[0] - 1)**2 + (x[1] - 2.5)**2

cons = (
    {"type": "ineq", "fun": lambda x:  x[0] - 2*x[1] + 2},
    {"type": "ineq", "fun": lambda x: -x[0] - 2*x[1] + 6},
    {"type": "ineq", "fun": lambda x: -x[0] + 2*x[1] + 2}
)

bnds = ( (0, None), (0, None) )

x0 = [2.0, 0.0]
result = minimize(fun, x0, method="SLSQP", bounds=bnds, constraints=cons)

print("result = ")
print(result)
