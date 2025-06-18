using GLMakie  # Using GLMakie for 3D visualization
using Colors
using LinearAlgebra
using Statistics: mean

"""
    visualize_mesh(points, faces; boundary_data=nothing)

Visualize an OpenFOAM mesh with optional boundary coloring.
"""
function visualize_mesh(points, faces; boundary_data=nothing, owner=nothing, boundary=nothing)
    fig = Figure(resolution=(1200, 800))
    
    # Create 3D axis
    ax = Axis3(fig[1, 1], 
              title="OpenFOAM Mesh Visualization",
              perspectiveness=0.5,
              azimuth=0.5π,
              elevation=0.3π)
    
    # Calculate mesh center and scale for proper viewing
    mesh_center = mean(points, dims=1)
    mesh_size = maximum([norm(p .- mesh_center) for p in eachrow(points)])
    
    # Create colormap for boundaries
    boundary_colors = if boundary_data !== nothing
        distinguishable_colors(length(boundary_data), [RGB(1,1,1)])
    else
        [RGBA(0.8, 0.8, 0.8, 0.6)] # Default color
    end
    
    # Plot all faces
    for (i, face) in enumerate(faces)
        face_points = points[face, :]
        
        # Determine color (check if face is a boundary face)
        color = if boundary_data !== nothing && owner !== nothing && boundary !== nothing
            # Check if this face is a boundary face
            boundary_idx = findfirst(b -> owner[i] in b["faceRange"], values(boundary))
            boundary_idx !== nothing ? boundary_colors[boundary_idx] : RGBA(0.8, 0.8, 0.8, 0.2)
        else
            RGBA(0.8, 0.8, 0.8, 0.6)
        end
        
        # Plot the face polygon
        poly!(ax, face_points, color=color, strokecolor=:black, strokewidth=0.5, transparency=true)
        
        # Add wireframe
        lines!(ax, [face_points; face_points[1:1,:]], color=:black, linewidth=0.5)
    end
    
    # Add boundary legend if available
    if boundary_data !== nothing
        leg_labels = [Text("$(k)") for k in keys(boundary_data)]
        leg_colors = boundary_colors[1:length(boundary_data)]
        
        Legend(fig[1, 2], 
               [PolyElement(color=c) for c in leg_colors],
               leg_labels,
               "Boundaries")
    end
    
    # Set axis limits based on mesh dimensions
    xlims!(ax, mesh_center[1]-mesh_size, mesh_center[1]+mesh_size)
    ylims!(ax, mesh_center[2]-mesh_size, mesh_center[2]+mesh_size)
    zlims!(ax, mesh_center[3]-mesh_size, mesh_center[3]+mesh_size)
    
    # Add helpful controls
    ax.scene.camera.fov = 30  # Field of view
    display(fig)
    
    return fig
end

"""
    visualize_boundaries(points, faces, owner, boundary)

Specialized visualization showing only boundary faces with coloring.
"""
function visualize_boundaries(points, faces, owner, boundary)
    # Create boundary face ranges
    boundary_face_ranges = Dict{String, UnitRange{Int}}()
    for (name, data) in boundary
        start = data["startFace"] + 1  # 1-based indexing
        stop = start + data["nFaces"] - 1
        boundary_face_ranges[name] = start:stop
    end
    
    fig = Figure(resolution=(1200, 800))
    ax = Axis3(fig[1, 1], title="OpenFOAM Boundary Faces")
    
    # Create colormap for boundaries
    boundary_colors = distinguishable_colors(length(boundary), [RGB(1,1,1)])
    
    # Plot boundary faces
    for (i, (name, range)) in enumerate(boundary_face_ranges)
        for face_idx in range
            face = faces[face_idx]
            face_points = points[face, :]
            
            poly!(ax, face_points, 
                 color=(boundary_colors[i], 0.7),
                 strokecolor=:black,
                 strokewidth=0.5)
            
            #lines!(ax, [face_points; face_points[1:1,:]], 
            #      color=:black, 
            #      linewidth=0.5)
        end
    end
    
    # Add legend
    Legend(fig[1, 2], 
           [PolyElement(color=c) for c in boundary_colors],
           collect(keys(boundary)),
           "Boundaries")
    
    # Set nice view
    mesh_center = mean(points, dims=1)
    mesh_size = maximum([norm(p .- mesh_center) for p in eachrow(points)])
    
    xlims!(ax, mesh_center[1]-mesh_size, mesh_center[1]+mesh_size)
    ylims!(ax, mesh_center[2]-mesh_size, mesh_center[2]+mesh_size)
    zlims!(ax, mesh_center[3]-mesh_size, mesh_center[3]+mesh_size)
    
    display(fig)
    return fig
end

"""
    simple_wireframe(points, faces)

Minimal wireframe visualization for quick viewing.
"""
function simple_wireframe(points, faces)
    fig = Figure()
    ax = Axis3(fig[1, 1], title="Mesh Wireframe")
    
    for face in faces
        face_points = points[face, :]
        lines!(ax, [face_points; face_points[1:1,:]], 
              color=(:black, 0.5),
              linewidth=0.3)
    end
    
    # Set equal aspect ratio
    mesh_center = mean(points, dims=1)
    mesh_size = maximum([norm(p .- mesh_center) for p in eachrow(points)])
    
    xlims!(ax, mesh_center[1]-mesh_size, mesh_center[1]+mesh_size)
    ylims!(ax, mesh_center[2]-mesh_size, mesh_center[2]+mesh_size)
    zlims!(ax, mesh_center[3]-mesh_size, mesh_center[3]+mesh_size)
    
    display(fig)
    return fig
end

