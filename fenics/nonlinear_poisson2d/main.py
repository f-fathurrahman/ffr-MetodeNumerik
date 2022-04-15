from fenics import *

def q(u):
    return 1 + u**2

u_code = "x[0] + 2*x[1] + 1"
f_code = "-10*x[0] - 20*x[1] - 10"

mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, "P", 1)
u_D = Expression(u_code, degree=1)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

u = Function(V)
v = TestFunction(V)
f = Expression(f_code, degree=1)

F = q(u)*dot( grad(u), grad(v) )*dx - f*v*dx

solve(F == 0, u, bc)

vtkfile = File("OUT_solution.pvd")
vtkfile << u

