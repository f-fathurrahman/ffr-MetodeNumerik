using Revise
using Infiltrator
using DelimitedFiles

includet("INC_read_polyMesh_v2.jl")
includet("INC_plotMesh_v1.jl")

# Load your mesh data first (using the previous reading functions)
mesh_dir = "constant/polyMesh"
points = read_points(joinpath(mesh_dir, "points"))
faces = read_faces(joinpath(mesh_dir, "faces"))
owner = read_owner_neighbor(joinpath(mesh_dir, "owner"))
boundary = read_boundary(joinpath(mesh_dir, "boundary"))

# 1. Full mesh visualization with boundary coloring
#fig1 = visualize_mesh(points, faces, owner=owner, boundary=boundary)

# 2. Boundary-only visualization
fig2 = visualize_boundaries(points, faces, owner, boundary)

# 3. Quick wireframe view
#fig3 = simple_wireframe(points, faces)

