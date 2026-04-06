import meshio
import numpy as np
from my_structs import MagnetostaticElement
from my_solvers import MagnetostaticSolver

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.patches import Polygon

# Load mesh and extract geometry
mesh = meshio.read("iron_core.msh")
points = mesh.points[:,:2] # Extract XY coordinates
triangles = mesh.cells_dict["triangle"] # Triangular elements
tags = mesh.cell_data_dict["gmsh:physical"]["triangle"] # Tags for triangles
line_tags = mesh.cell_data_dict["gmsh:physical"]["line"] # Tags for lines
lines = mesh.cells_dict["line"] # Line elements
# Map physical tags to names
map_f = {
    v[0]: k for k, v in mesh.field_data.items()
}

# Material properties
materials = {
    "Iron": {
        "nu": 1 / 1000.0,
        "J": 0.0
    }, # Magnetic permeability and current density
    "Coil+": {
        "nu": 1,
        "J": 1e4
    }, # Positive current density for the coil
    "Coil-": {
        "nu": 1,
        "J": -1e4
    } # Negative current density for the coil
}


# Initialize elements and solver
elements = []
nodeTags = np.arange(len(points)) # Node indices
triElements = []
# Loop to create magnetostatic elements
for i, tri in enumerate(triangles):
    name = map_f.get(tags[i]) # Get the region name
    prop = materials.get(name, {"nu": 1, "J": 0.0}) # Default properties if not found
    x, y = points[tri, 0], points[tri, 1] # Node coordinates
    # Instantiate the magnetostatic element
    e = MagnetostaticElement()
    e.setNodes(x, y)
    e.setProperties(prop["nu"], prop["J"])
    elements.append(e)
    triElements.append(tri)

# Define boundary conditions
boundary_conditions = {
    "A0": {
        "nodes": [],
        "potential": 0.0
    } # Zero potential on boundary "A0"
}

# Map lines to physical groups
bc_tag = mesh.field_data["A0"][0] # Tag for boundary "A0"
for line, tag in zip(lines, line_tags):
    if tag == bc_tag:
        boundary_conditions["A0"]["nodes"].extend(line)

# Initialize the solver
solver = MagnetostaticSolver(
    nodes=points,
    triElements=triElements,
    elements=elements,
    boundaryConditions=boundary_conditions
)

# Assemble the global matrix and force vector
solver.assemble_global_matrix_and_vector()

# Apply boundary conditions
solver.apply_boundary_conditions()

# Solve the linear system
solver.solve()

# Calculate the magnetic field
solver.calculate_magnetic_field()

# Get the results
V = solver.get_potential() # Magnetic potential
Bx, By, modB, xc, yc = solver.get_magnetic_field() # Magnetic field components and centroids

# Visualization
triang = mtri.Triangulation(points[:, 0], points[:, 1], triangles) # Create triangulation
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
# Iron polygon points (from .geo file)
iron_points = [
    (-17, -24.8), (23, -24.8), (23, -0.3), (13, -0.3),
    (13, -15.3), (-7, -15.3), (-7, 14.7), (13, 14.7),
    (13, 1.7), (23, 1.7), (23, 24.7), (-17, 24.7)
]

# Plot Magnetic Potential A
pot = axes[0].tricontourf(triang, np.abs(V), levels=50, cmap="coolwarm")
axes[0].tricontour(triang, np.abs(V), levels=25, colors="black", linewidths=0.5)
fig.colorbar(pot, ax=axes[0], fraction=0.1, pad=0.04, label="Potential A")
axes[0].set_title("Magnetic Potential A")
axes[0].set_xlabel("x (m)")
axes[0].set_ylabel("y (m)")
axes[0].set_aspect("equal")

# Add iron core polygon outline
iron_poly2 = Polygon(iron_points, closed=True, edgecolor="k", facecolor="none", linewidth=1)
axes[0].add_patch(iron_poly2)

# Plot Magnetic Field
field = axes[1].tripcolor(triang, modB, shading="flat", cmap="coolwarm")
fig.colorbar(field, ax=axes[1], fraction=0.1, pad=0.04, label="|B| (T)")

# Uncomment for vector field
# axes[1].quiver(xc, yc, Bx, By, scale=20, color="white", width=0.002, alpha=0.8)
axes[1].set_title("Magnetic Field Magnitude and Direction")
axes[1].set_xlabel("x (m)")
axes[1].set_ylabel("y (m)")
axes[1].set_aspect("equal")

# Add iron core polygon outline
iron_poly2 = Polygon(iron_points, closed=True, edgecolor="k", facecolor="none", linewidth=1)
axes[1].add_patch(iron_poly2)
plt.tight_layout()
plt.savefig("IMG_magn_field.png", dpi=300, bbox_inches="tight")
plt.show()

