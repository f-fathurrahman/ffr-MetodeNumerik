from fenics import *
from ufl import nabla_div

L = 1
W = 0.2
mu = 1
rho = 1
delta = W/L
gamma = 0.4*delta**2
beta = 1.25
lambda_ = beta
g = gamma

# Mesh
mesh = BoxMesh(Point(0,0,0), Point(L,W,W), 10, 3, 3)
V = VectorFunctionSpace(mesh, "P", 1)

# Boundary condition
TOL = 1e-14

def clamped_boundary(x, on_boundary):
    return on_boundary and x[0] < TOL

bc = DirichletBC( V, Constant( (0,0,0) ), clamped_boundary )

# Strain and stress
def epsilon(u):
    return 0.5*( nabla_grad(u) + nabla_grad(u).T )


def sigma(u):
    return lambda_ * nabla_div(u)*Identity(d) + 2*mu*epsilon(u)

# Variational problem
u = TrialFunction(V)
d = u.geometric_dimension()
v = TestFunction(V)
f = Constant( (0,0,-rho*g) )
T = Constant( (0,0,0) )
a = inner( sigma(u), epsilon(v) )*dx
L = dot(f, v)*dx + dot(T, v)*ds

# Compute solution
u = Function(V)

solve( a == L, u, bc )

vtkfile = File("OUT_u.pvd")
vtkfile << u

# Calculate stress
s = sigma(u) - (1.0/3.0)*tr( sigma(u) ) * Identity(d)
von_Mises = sqrt( 3.0/2.0*inner(s,s) )
V = FunctionSpace(mesh, "P", 1)
von_Mises = project(von_Mises, V)

vtkfile = File("OUT_von_Mises_stress.pvd")
vtkfile << von_Mises

u_magnitude = sqrt( dot(u, u) )
u_magnitude = project(u_magnitude, V)

vtkfile = File("OUT_u_magnitude.pvd")
vtkfile << u_magnitude
