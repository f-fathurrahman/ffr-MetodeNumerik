from sympy import *

h = Symbol("h", real=True)
p = Symbol("p", real=True)

# N: number of interior nodes
def generate_trid_matrix(N):
    assert(N >= 1)
    A = Matrix.zeros(N)
    for i in range(N):
        A[i,i] = 2 - h**2 * p**2
        if i >= 1:
            A[i-1,i] = -1
            A[i,i-1] = -1
    return A

for i in range(1,5):
    print()
    A = generate_trid_matrix(i)
    A_subsh = A.subs({h: 3.0/(i+1)})
    #pprint(A)
    determinant = A_subsh.det()
    pprint(determinant)
    sols = solve(determinant)
    pprint(sols)
