import meshio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from my_structs import ElectrostaticElement
from my_solvers import ElectrostaticSolver

# 1. Load the mesh generated in Gmsh
data = meshio.read("square2.msh")
points = data.points[:,:2] # Node coordinates
triangles = data.cells_dict["triangle"] # Triangular elements
lines = data.cells_dict["line"] # Lines for boundary conditions

# 2. Map physical tags (from Gmsh) to names
physical_map = {v[0]: k for k, v in data.field_data.items()}
line_tags = data.cell_data_dict["gmsh:physical"]["line"]

# 3. Initialize electrostatic elements
elements = [ElectrostaticElement() for _ in range(len(triangles))]
for element, tri in zip(elements, triangles):
    element.setNodes(*points[tri].T)

# 4. Define boundary conditions
boundary_conditions = {
    "Top": {
        "nodes": [],
        "potential": 100.0
    },
    "Bottom": {
        "nodes": [],
        "potential": -100.0
    },
    "Left": {
        "nodes": [],
        "potential": 0.0
    },
    "Right": {
        "nodes": [],
        "potential": 0.0
    }
}

# 5. Associate boundary nodes based on Gmsh data
for line, tag in zip(lines, line_tags):
    if physical_map.get(tag) in boundary_conditions:
        boundary_conditions[physical_map[tag]]["nodes"].extend(line)

# 6. Remove duplicate nodes from boundary conditions
for bc in boundary_conditions.values():
    bc["nodes"] = list(set(bc["nodes"]))

# 7. Initialize the solver and solve the electrostatic system
solver = ElectrostaticSolver(
    points,
    triangles,
    elements,
    boundary_conditions
)
solver.assemble_global_matrix_and_vector()
solver.apply_boundary_conditions()
solver.solve()

# 8. Get the potential results
V = solver.get_potential()

# 9. Calculate the electric field
#solver.calculate_electric_field()
solver.solve()
Ex, Ey, modE, xc, yc = solver.get_electric_field()

# 10. Create the triangulation
triang = mtri.Triangulation(points[:, 0], points[:, 1], triangles)

# 11. Obtain the electric fields and centroids" coordinates
Ex, Ey, modE, xc, yc = solver.get_electric_field()


# 12. Prepare the subplots
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
sampling_factor = 1 # Denser
scale_factor = np.max(modE) * 3 # Increased for clear visibility
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
