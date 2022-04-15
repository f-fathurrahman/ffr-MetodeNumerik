from fenics import *
from mshr import *
import matplotlib.pyplot as plt

import matplotlib.style
matplotlib.style.use("classic")

# Mesh
domain = Circle(Point(0,0), 1)
mesh = generate_mesh(domain, 64)

# Load
beta = 8
R0 = 0.6
p = Expression("4*exp( -pow(beta,2)*( pow(x[0],2) + pow( x[1] - R0, 2) ) )", degree=1, beta=beta, R0=R0)

# Function space
V = FunctionSpace(mesh, "P", 1)

# Boundary conditione
w_D = Constant(0)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, w_D, boundary)

# Variational problem
w = TrialFunction(V)
v = TestFunction(V)
a = dot( grad(w), grad(v) )*dx
L = p*v*dx

# Solve the variational problem
w = Function(V)
solve(a == L, w, bc)

# Plot
p_V = interpolate(p, V) # pressure/Load, must be interpolated in this function space

plot(p_V, title="Load")
plt.savefig("IMG_load.png", dpi=150)

plot(w, title="Deflection")
plt.savefig("IMG_deflection.png", dpi=150)
