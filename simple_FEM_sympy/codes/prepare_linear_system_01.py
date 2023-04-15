exec(open("setup_problem_01.py").read())

# Functions, symbolic
T = Function("T")(x)
w = Function("w")(x)
k = Symbol("k", real=True, positive=True)
Q = Symbol("Q", real=True)


Tnodal = []
for i in range(Nnodes):
    Tnodal.append(Symbol("T_" + str(i), real=True))

T_expansion = 0
for i in range(Nnodes):
    T_expansion += Tnodal[i]*Nfuncs[i]

# Some concrete values
pprint(T_expansion.subs({x: 0.1, L: Lnum}))


# Stiffness term
expr1 = k * Derivative(w,x) * Derivative(T, x)
expr1s = expr1.subs({T: T_expansion, w: Nfuncs[0]})
pprint(expr1s)

expr1ss = expr1s.doit()
pprint(expr1ss)

I1 = Integral(expr1ss, (x,0,L))
I1

term1_first = I1.doit()


# Source term
expr2 = Integral( Q*w, (x,0,L) )
term2_first = expr2.subs({w: Nfuncs[0]}).doit()
term2_first


expr3 = k*w*Derivative(T,x)
expr3

term3_x0 = expr3.subs({ w: Nfuncs[0], x: 0, L: Lnum})
term3_x0


term3_xL = expr3.subs({ w: Nfuncs[0], x: Lnum, L: Lnum})
term3_xL

term3_first = term3_xL - term3_x0
term3_first

#

eq_first = term1_first - term2_first - term3_first
eq_first
q0 = Symbol("q_0", real=True)
term3_first
term3_first.args
term3_first.args[2]
#eq_first = eq_first.subs( {k*term3_first.args[2]: -q0} )
#eq_first
