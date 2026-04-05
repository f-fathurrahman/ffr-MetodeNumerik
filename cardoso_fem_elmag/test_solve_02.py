import meshio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

from my_structs import ElectrostaticElement
from my_solvers import ElectrostaticSolver

# Load mesh and extract geometry
mesh = meshio.read("parallel_plate.msh")
points = mesh.points[:, :2]
triangles = mesh.cells_dict["triangle"]
lines = mesh.cells_dict["line"]

# Physical group mapping
physical_map = {
    v[0]: k for k, v in mesh.field_data.items()
}
triangle_tags = mesh.cell_data_dict["gmsh:physical"]["triangle"]
line_tags = mesh.cell_data_dict["gmsh:physical"]["line"]

# Material assignment based on physical regions
elements = []
nodeTags = np.arange(len(points))
triElements = []
# Permittivity by region
region_eps = {
    "Insulator": 3.5,
    "External": 1.0
}

# Loop through elements and assign material
for tri, tag in zip(triangles, triangle_tags):
    name = physical_map.get(tag)
    eps = region_eps.get(name, 1.0)
    coords = points[tri]
    # Using the correct class
    e = ElectrostaticElement()
    e.setNodes(coords[:, 0], coords[:, 1])
    e.setProperties(eps=eps, rho=0.0)
    elements.append(e)
    triElements.append(tri)

# Dirichlet boundary conditions
boundary_conditions = {
    "V0": {
        "nodes": [],
        "potential": 2.0
    },
    "GND": {
        "nodes": [],
        "potential": -2.0
    }
}

# Mapping lines to Physical Groups
for line, tag in zip(lines, line_tags):
    name = physical_map.get(tag)
    if name in boundary_conditions:
        boundary_conditions[name]["nodes"].extend(line)

# Solver Initialization
solver = ElectrostaticSolver(
    nodes=points,
    triElements=triElements,
    elements=elements,
    boundaryConditions=boundary_conditions
)

# Assemble the global matrix and the force vector
solver.assemble_global_matrix_and_vector()
# Apply boundary conditions
solver.apply_boundary_conditions()
# Solve the linear system
solver.solve()
# Calculate the electric fields
#solver.calculate_electric_field()

# Obtain the results
V = solver.get_potential()
Ex, Ey, modE, xc, yc = solver.get_electric_field()

triang = mtri.Triangulation(points[:, 0], points[:, 1], triangles)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# ---------- Subplot 1 | Electrostatic Potential ----------
pot = axes[0].tricontourf(triang, V, levels=50, cmap="jet")
fig.colorbar(pot, ax=axes[0], label="Potential (V)")
axes[0].set_title("Electrostatic Potential")
axes[0].set_xlabel("x (m)")
axes[0].set_ylabel("y (m)")
axes[0].set_aspect("equal")

# ---------- Subplot 2 | Electric Field ----------
# Smooth interpolation of the field
field = axes[1].tripcolor(triang, modE, shading="flat", cmap="jet")
fig.colorbar(field, ax=axes[1], label="|E| (V/m)")
# Adjust scale and density
sampling_factor = 2 # smaller -> Denser
scale_factor = np.max(modE) * 30 # Increased for clear visibility
axes[1].quiver(xc[::sampling_factor], yc[::sampling_factor],
               Ex[::sampling_factor], Ey[::sampling_factor],
               scale=scale_factor, color="white", width=0.0025)
axes[1].set_title("Electric Field Magnitude and Vectors")
axes[1].set_xlabel("x (m)")
axes[1].set_ylabel("y (m)")
axes[1].set_aspect("equal")

# ---------- Finalization ----------
plt.tight_layout()
plt.show()

