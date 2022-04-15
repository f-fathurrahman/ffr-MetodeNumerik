from fenics import *
import numpy as np
import matplotlib.pyplot as plt

import matplotlib.style
matplotlib.style.use("classic")

T = 2.0 # final time
Nsteps = 10 # time steps
dt = T/Nsteps
alpha = 3
beta = 1.2

# Mesh
Nx = Ny = 8
mesh = UnitSquareMesh(Nx, Ny)

# Function space
V = FunctionSpace(mesh, "P", 1)

# Boundary condition
u_D = Expression("1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t",
    degree=2, alpha=alpha, beta=beta, t=0)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Initial value
u_n = interpolate(u_D, V)

# Variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(beta - 2 - 2*alpha)

F = u*v*dx + dt*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx
a, L = lhs(F), rhs(F)

# Time stepping
u = Function(V)
t = 0.0
vtkfile = File("output/solution.pvd")

for n in range(Nsteps):
    #
    t = t + dt
    u_D.t = t
    # Compute solution
    solve(a == L, u, bc)
    #
    
    #plot(u)
    #plt.savefig("IMG_u_{:03d}".format(n+1))
    vtkfile << (u, t)

    # Compute error at vertices
    u_e = interpolate(u_D, V)
    #
    vertex_values_u_e = u_e.compute_vertex_values(mesh)
    vertex_values_u = u.compute_vertex_values(mesh)

    error = np.max( np.abs(vertex_values_u_e - vertex_values_u) )
    print('t = %.2f: error = %.3g' % (t, error))
    # Update previous solution
    u_n.assign(u)

