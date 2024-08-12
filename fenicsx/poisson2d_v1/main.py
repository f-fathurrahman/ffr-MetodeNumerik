import dolfinx

import petsc4py.PETSc
from mpi4py import MPI

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




