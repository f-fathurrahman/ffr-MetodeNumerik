using Infiltrator
using DelimitedFiles


"""
    read_points(filepath)

Read OpenFOAM points file and return coordinates as Matrix{Float64}.
"""
function read_points(filepath)
    # Read the entire file
    content = read(filepath, String)
    
    # Find the content between the parentheses
    start_idx = findfirst('(', content)
    end_idx = findlast(')', content)
    points_content = content[start_idx+1:end_idx-1]
    
    # Parse the numbers
    points = []
    for line in split(points_content, '\n')
        # Remove comments and trim whitespace
        clean_line = split(line, "//")[1]
        clean_line = strip(clean_line)
        
        if !isempty(clean_line)
            # Remove parentheses and split numbers
            coords = split(replace(clean_line, r"[()]" => ""))
            push!(points, parse.(Float64, coords))
        end
    end
    
    # Convert to matrix (N×3)
    return vcat(points'...)
end


"""
    read_faces(filepath)

Read OpenFOAM faces file in compact format (e.g., "4(0 1 2 3)") 
and return as Vector{Vector{Int64}}.
"""
function read_faces(filepath)
    content = read(filepath, String)
    
    # Find the content between the main parentheses
    start_idx = findfirst('(', content)
    end_idx = findlast(')', content)
    faces_content = content[start_idx+1:end_idx-1]
    
    faces = Vector{Vector{Int64}}()
    current_face = Int64[]
    
    # Use a more robust parsing approach for compact format
    for m in eachmatch(r"(\d+)\(([\d\s]+)\)", faces_content)
        n_points = parse(Int64, m.captures[1])
        points_str = m.captures[2]
        
        # Parse the vertex indices
        vertices = parse.(Int64, split(points_str)) .+ 1  # Convert to 1-based
        
        # Verify we got the correct number of points
        if length(vertices) != n_points
            @warn "Face point count mismatch: expected $n_points, got $(length(vertices))"
        end
        
        push!(faces, vertices)
    end
    
    # Fallback if regex fails (shouldn't happen with proper OpenFOAM files)
    if isempty(faces)
        @warn "Compact format parsing failed, trying alternative approach"
        return read_faces_fallback(faces_content)
    end
    
    return faces
end

"""
    read_faces_fallback(content)

Fallback face reader for non-standard formats.
"""
function read_faces_fallback(content)
    faces = Vector{Vector{Int64}}()
    current_face = Int64[]
    expecting_count = true
    n_points = 0
    
    for token in split(content)
        if expecting_count
            # First number is the point count
            n_points = parse(Int64, token)
            expecting_count = false
        else
            # Subsequent numbers are vertex indices (might have trailing ')')
            num = parse(Int64, replace(token, r"[^0-9]" => ""))
            push!(current_face, num + 1)  # 0-based to 1-based
            
            if length(current_face) == n_points
                push!(faces, current_face)
                current_face = Int64[]
                expecting_count = true
            end
        end
    end
    
    return faces
end


"""
    read_owner_neighbor(filepath)

Read owner or neighbor file and return as Vector{Int64}.
"""
function read_owner_neighbor(filepath)
    content = read(filepath, String)
    
    # Find the content between the parentheses
    start_idx = findfirst('(', content)
    end_idx = findlast(')', content)
    data_content = content[start_idx+1:end_idx-1]
    
    data = Int64[]
    
    for line in split(data_content, '\n')
        # Remove comments and trim whitespace
        clean_line = split(line, "//")[1]
        clean_line = strip(clean_line)
        
        if !isempty(clean_line)
            push!(data, parse(Int64, clean_line) + 1)  # Convert to 1-based indexing
        end
    end
    
    return data
end

"""
    read_boundary(filepath)

Read boundary file and return as Dict with boundary information.
"""
function read_boundary(filepath)
    content = read(filepath, String)
    
    # Find the content between the parentheses
    start_idx = findfirst('(', content)
    end_idx = findlast(')', content)
    boundary_content = content[start_idx+1:end_idx-1]
    
    boundaries = Dict()
    current_key = ""
    current_dict = Dict()
    in_block = false
    brace_count = 0
    
    for line in split(boundary_content, '\n')
        # Remove comments and trim whitespace
        clean_line = split(line, "//")[1]
        clean_line = strip(clean_line)
        
        if isempty(clean_line)
            continue
        end
        
        if !in_block
            # Looking for boundary name
            if !occursin(r"[{}]", clean_line)
                current_key = clean_line
                in_block = true
                brace_count = 0
                current_dict = Dict()
            end
        else
            # Inside a boundary block
            if occursin('{', clean_line)
                brace_count += count(c -> c == '{', clean_line)
            end
            
            if occursin('}', clean_line)
                brace_count -= count(c -> c == '}', clean_line)
            end
            
            # Parse key-value pairs
            if occursin(';', clean_line)
                parts = split(replace(clean_line, ";" => ""))
                key = parts[1]
                value = join(parts[2:end], " ")
                
                # Try to parse as Int or Float if possible
                try
                    value = parse(Int64, value)
                catch
                    try
                        value = parse(Float64, value)
                    catch
                        # Keep as string
                    end
                end
                
                current_dict[key] = value
            end
            
            # End of block
            if brace_count == 0
                boundaries[current_key] = current_dict
                in_block = false
            end
        end
    end
    
    return boundaries
end


function debug_main()
    # Example usage:
    mesh_dir = joinpath("constant", "polyMesh")

    points = read_points(joinpath(mesh_dir, "points"))
    faces = read_faces(joinpath(mesh_dir, "faces"))
    owner = read_owner_neighbor(joinpath(mesh_dir, "owner"))
    neighbor = read_owner_neighbor(joinpath(mesh_dir, "neighbour"))
    boundary = read_boundary(joinpath(mesh_dir, "boundary"))

    println("Read mesh with $(size(points, 1)) points, $(length(faces)) faces, and $(length(boundary)) boundaries.")

    @infiltrate

end


