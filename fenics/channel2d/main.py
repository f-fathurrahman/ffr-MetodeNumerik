from fenics import *
import numpy as np

T = 10.0
Nsteps = 500
dt = T/Nsteps

mu = 1 # kinematic viscosity
rho = 1 # density

# Mesh
mesh = UnitSquareMesh(16, 16)

# Function spaces
V = VectorFunctionSpace(mesh, "P", 2)
Q = FunctionSpace(mesh, "P", 1)

# Boundaries
inflow = "near(x[0], 0)"
outflow = "near(x[0], 1)"
walls = "near(x[1], 0) || near(x[1], 1)"

# Boundary conditions
bcu_noslip = DirichletBC(V, Constant((0,0)), walls)
bcp_inflow = DirichletBC(Q, Constant(8), inflow)
bcp_outflow = DirichletBC(Q, Constant(0), outflow)
bcu = [bcu_noslip]
bcp = [bcp_inflow, bcp_outflow]

# Trial and test functions
u = TrialFunction(V)
v = TestFunction(V)

p = TrialFunction(Q)
q = TestFunction(Q)

# Functions for solutions at previous and current time steps
u_n = Function(V)
u_ = Function(V)

p_n = Function(Q)
p_ = Function(Q)

# Expressions used in variational forms
U = 0.5*(u_n + u)
n = FacetNormal(mesh)
f = Constant((0,0))
k = Constant(dt)
mu = Constant(mu)
rho = Constant(rho)

# Strain-rate tensor
def epsilon(u):
    return sym(nabla_grad(u))

def sigma(u, p):
    return 2*mu*epsilon(u) - p*Identity(len(u))

# Define variational problem for step 1
F1 = rho*dot((u - u_n) / k, v)*dx + \
     rho*dot(dot(u_n, nabla_grad(u_n)), v)*dx \
   + inner(sigma(U, p_n), epsilon(v))*dx \
   + dot(p_n*n, v)*ds - dot(mu*nabla_grad(U)*n, v)*ds \
   - dot(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Define variational problem for step 2
a2 = dot(nabla_grad(p), nabla_grad(q))*dx
L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx - (1/k)*div(u_)*q*dx

# Define variational problem for step 3
a3 = dot(u, v)*dx
L3 = dot(u_, v)*dx - k*dot(nabla_grad(p_ - p_n), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Apply boundary conditions to matrices
[bc.apply(A1) for bc in bcu]
[bc.apply(A2) for bc in bcp]

#vtkfile_u = File("output/OUT_u.pvd")

t = 0
for n in range(Nsteps):

    # Update current time
    t = t + dt

    # Step 1: Tentative velocity step
    b1 = assemble(L1)
    [bc.apply(b1) for bc in bcu]
    solve(A1, u_.vector(), b1)

    # Step 2: Pressure correction step
    b2 = assemble(L2)
    [bc.apply(b2) for bc in bcp]
    solve(A2, p_.vector(), b2)

    # Step 3: Velocity correction step
    b3 = assemble(L3)
    solve(A3, u_.vector(), b3)

    # Plot solution
    #vtkfile_u << u_

    # Compute error
    u_e = Expression( ("4*x[1]*(1.0 - x[1])", "0"), degree=2 )
    u_e = interpolate(u_e, V)
    error = np.abs( u_e.vector() - u_.vector() ).max()
    print("t = %18.10f: error = %18.10f" % (t, error))
    print("max u:", u_.vector().max())

    # Update previous solution
    u_n.assign(u_)
    p_n.assign(p_)

