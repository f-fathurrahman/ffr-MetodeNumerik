using Revise
using Infiltrator
using DelimitedFiles

includet("INC_read_polyMesh_v2.jl")


function debug_main()
    mesh_dir = "constant/polyMesh"

    # Read with optimized functions
    points = read_points(joinpath(mesh_dir, "points"))  # N×3 Matrix
    faces = read_faces(joinpath(mesh_dir, "faces"))     # Vector{Vector{Int}}
    owner = read_owner_neighbor(joinpath(mesh_dir, "owner"))
    neighbor = read_owner_neighbor(joinpath(mesh_dir, "neighbour"))
    boundary = read_boundary(joinpath(mesh_dir, "boundary"))

    println("Mesh statistics:")
    println("- Points: ", size(points, 1))
    println("- Faces: ", length(faces))
    println("- Owner cells: ", length(owner))
    println("- Neighbor cells: ", length(neighbor))
    println("- Boundaries: ", length(boundary))

    @exfiltrate
end


