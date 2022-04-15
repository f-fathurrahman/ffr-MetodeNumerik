from fenics import *
import numpy as np
import matplotlib.pyplot as plt

T = 2.0 # final time
Nsteps = 50 # time steps
dt = T/Nsteps

# Mesh
Nx = Ny = 30
mesh = RectangleMesh(Point(-2,-2), Point(2,2), Nx, Ny)

# Function space
V = FunctionSpace(mesh, "P", 1)

# Boundary condition
def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, Constant(0.0), boundary)

# Initial value
u_0 = Expression("exp(-a*pow(x[0],2) - a*pow(x[1],2))", degree=2, a=5)
# interpolate to set initial value
u_n = interpolate(u_0, V)

# Variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0.0)

F = u*v*dx + dt*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx
a, L = lhs(F), rhs(F)

# Time stepping
u = Function(V)
t = 0.0
vtkfile = File("output/solution.pvd")

for n in range(Nsteps):

    t = t + dt

    # Compute solution
    solve(a == L, u, bc)

    vtkfile << (u, t)
    #
    # Update previous solution
    u_n.assign(u)

