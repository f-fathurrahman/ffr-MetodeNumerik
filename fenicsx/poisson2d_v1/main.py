# From:
# https://docs.fenicsproject.org/dolfinx/v0.8.0/python/demos/demo_poisson.html

import dolfinx

import petsc4py.PETSc
from mpi4py import MPI
from petsc4py.PETSc import ScalarType

import numpy as np
import ufl
from dolfinx import fem, io, mesh, plot
from dolfinx.fem.petsc import LinearProblem
from ufl import ds, dx, grad, inner

# Create a rectangular mesh
msh = mesh.create_rectangle(
    comm=MPI.COMM_WORLD,
    points=( (0.0, 0.0), (2.0, 1.0) ),
    n=(32,16),
    cell_type=mesh.CellType.triangle
)

# Function space
V = fem.functionspace(msh, ("Lagrange", 1))

# Find the boundary
facets = mesh.locate_entities_boundary(
    msh,
    dim=(msh.topology.dim - 1),
    marker=lambda x: np.isclose(x[0], 0.0) | np.isclose(x[0], 2.0),
)

# Find the DOS that are associated with boundary facets
dofs = fem.locate_dofs_topological(V=V, entity_dim=1, entities=facets)
bc = fem.dirichletbc(value=ScalarType(0), dofs=dofs, V=V)

# Define variational problem
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
x = ufl.SpatialCoordinate(msh)

f = 10 * ufl.exp(-((x[0] - 0.5) ** 2 + (x[1] - 0.5) ** 2) / 0.02) # RHS function
g = ufl.sin(5*x[0])
a = inner(grad(u), grad(v)) * dx
L = inner(f, v) * dx + inner(g, v)*ds

# Create LinearProblem from the variational problem and BC
problem = LinearProblem(
    a, L, bcs=[bc],
    petsc_options={
        "ksp_type": "preonly",
        "pc_type": "lu"
    }
)
uh = problem.solve()

# Save solution to XDM file
with io.XDMFFile(msh.comm, "TEMP_poisson.xdmf", "w") as file:
    file.write_mesh(msh)
    file.write_function(uh)

# Using pyvista
import pyvista
cells, types, x = plot.vtk_mesh(V) # read mesh
grid = pyvista.UnstructuredGrid(cells, types, x)
grid.point_data["u"] = uh.x.array.real
grid.set_active_scalars("u")
plotter = pyvista.Plotter()
plotter.add_mesh(grid, show_edges=True)
warped = grid.warp_by_scalar()
plotter.add_mesh(warped)
plotter.show()

# if pyvista.OFF_SCREEN
#pyvista.start_xvfb(wait=0.1) # need to install something
#plotter.screenshot("IMG_uh_poisson.png")

