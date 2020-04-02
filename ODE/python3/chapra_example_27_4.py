from sympy import *

k = Symbol("k", real=True, positive=True)
ω = Symbol("omega", real=True, positive=True)
m_1 = Symbol("m_1", real=True, positive=True)
m_2 = Symbol("m_2", real=True, positive=True)

# Matrix elements
M = Matrix.zeros(2,2)
M[0,0] = 2*k/m_1 - ω**2
M[0,1] = -k/m_1
M[1,0] = -k/m_2
M[1,1] = (2*k/m_2 - ω**2)

pprint(M)

# Substitute the numerical values
M_numeric = M.subs({m_1: 40, m_2: 40, k: 200})
# Calc determinant
solve( M_numeric.det() )