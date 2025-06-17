using Infiltrator
using DelimitedFiles


"""
    read_points(filepath)

Read OpenFOAM points file with preallocation.
Returns N×3 Matrix{Float64}.
"""
function read_points(filepath)
    content = read(filepath, String)
    
    # More robust header parsing
    header_lines = split(content, '\n')
    n_points = 0
    for line in header_lines
        if occursin(r"^\d+$", strip(line))  # Find the line with just a number
            n_points = parse(Int, strip(line))
            break
        end
    end
    n_points == 0 && error("Could not determine number of points from file header")
    
    points = Matrix{Float64}(undef, n_points, 3)
    
    # Find the content between parentheses
    data_start = findfirst('(', content)
    data_end = findlast(')', content)
    data_start === nothing && error("No data found between parentheses")
    points_content = content[data_start+1:data_end-1]
    
    # Parse all points at once
    point_entries = split(points_content, '\n', keepempty=false)
    point_idx = 1
    
    for entry in point_entries
        coords = match(r"\(([^)]+)\)", entry)
        if coords !== nothing
            points[point_idx,:] = parse.(Float64, split(coords.captures[1]))
            point_idx += 1
        elseif occursin(r"\d", entry)  # Fallback for formats without parentheses
            points[point_idx,:] = parse.(Float64, split(entry))
            point_idx += 1
        end
        point_idx > n_points && break
    end
    
    return points
end

"""
    read_faces(filepath)

Read OpenFOAM faces file with preallocation.
Returns Vector{Vector{Int}}.
"""
function read_faces(filepath)
    content = read(filepath, String)
    
    # Robust header parsing
    header_lines = split(content, '\n')
    n_faces = 0
    for line in header_lines
        if occursin(r"^\d+$", strip(line))
            n_faces = parse(Int, strip(line))
            break
        end
    end
    n_faces == 0 && error("Could not determine number of faces from file header")
    
    faces = Vector{Vector{Int}}(undef, n_faces)
    
    # Find the content between parentheses
    data_start = findfirst('(', content)
    data_end = findlast(')', content)
    data_start === nothing && error("No data found between parentheses")
    faces_content = content[data_start+1:data_end-1]
    
    # Process each face entry
    face_entries = split(faces_content, '\n', keepempty=false)
    face_idx = 1
    
    for entry in face_entries
        isempty(strip(entry)) && continue
        
        # Match compact format: "4(0 1 2 3)"
        m = match(r"(\d+)\(([\d\s]+)\)", entry)
        if m !== nothing
            n_points = parse(Int, m.captures[1])
            points = parse.(Int, split(m.captures[2])) .+ 1  # 0-based to 1-based
            faces[face_idx] = points
            face_idx += 1
        else
            # Fallback for non-compact format
            tokens = parse.(Int, split(entry))
            if !isempty(tokens)
                # Handle both formats: with/without leading count
                if length(tokens) > 1 && tokens[1] == length(tokens)-1
                    faces[face_idx] = tokens[2:end] .+ 1
                else
                    faces[face_idx] = tokens .+ 1
                end
                face_idx += 1
            end
        end
        face_idx > n_faces && break
    end
    
    return faces
end



"""
    read_owner_neighbor(filepath)

Read owner/neighbor files with preallocation.
Returns Vector{Int}.
"""
function read_owner_neighbor(filepath)
    content = read(filepath, String)
    
    # Robust header parsing
    header_lines = split(content, '\n')
    n_entries = 0
    for line in header_lines
        if occursin(r"^\d+$", strip(line))
            n_entries = parse(Int, strip(line))
            break
        end
    end
    n_entries == 0 && error("Could not determine number of entries from file header")
    
    data = Vector{Int}(undef, n_entries)
    
    # Find the content between parentheses
    data_start = findfirst('(', content)
    data_end = findlast(')', content)
    data_start === nothing && error("No data found between parentheses")
    entries = split(content[data_start+1:data_end-1], '\n', keepempty=false)
    
    idx = 1
    for entry in entries
        val = strip(split(entry, "//")[1])  # Remove comments
        if !isempty(val)
            data[idx] = parse(Int, val) + 1  # 0-based to 1-based
            idx += 1
        end
        idx > n_entries && break
    end
    
    return data
end



"""
    read_boundary(filepath)

Read OpenFOAM boundary file with proper dictionary parsing.
Returns Dict{String, Dict{String, Any}}.
"""
function read_boundary(filepath)
    content = read(filepath, String)
    
    # Extract number of boundaries from content
    n_boundaries = 0
    for line in split(content, '\n')
        stripped = strip(line)
        if occursin(r"^\d+$", stripped)
            n_boundaries = parse(Int, stripped)
            break
        end
    end
    
    # Find the main content between outer parentheses
    data_start = findfirst('(', content)
    data_end = findlast(')', content)
    (data_start === nothing || data_end === nothing) && error("No boundary data found between parentheses")
    boundary_content = strip(content[data_start+1:data_end-1])
    
    boundaries = Dict{String, Dict{String, Any}}()
    current_boundary = ""
    current_dict = Dict{String, Any}()
    brace_level = 0
    in_boundary = false
    
    # Process line by line
    for line in split(boundary_content, '\n')
        line = strip(split(line, "//")[1])  # Remove comments
        isempty(line) && continue
        
        if !in_boundary
            # Looking for boundary name
            if !isempty(line) && !occursin(r"[{};]", line)
                current_boundary = line
                in_boundary = true
                brace_level = 0
            end
        else
            # Track brace levels
            brace_level += count(c -> c == '{', line)
            
            # Parse key-value pairs
            if occursin(';', line)
                parts = split(replace(line, r"\s*;\s*" => ""))
                if length(parts) >= 2
                    key = parts[1]
                    value = join(parts[2:end], " ")
                    
                    # Try to parse as Int or Float if possible
                    try 
                        value = parse(Int, value)
                    catch
                        try 
                            value = parse(Float64, value)
                        catch
                            # Keep as string if not a number
                            # Remove quotes if present
                            value = replace(value, r"^\"|\"$" => "")
                        end
                    end
                    
                    current_dict[key] = value
                end
            end
            
            brace_level -= count(c -> c == '}', line)
            
            # End of boundary block
            if brace_level <= 0
                boundaries[current_boundary] = current_dict
                current_boundary = ""
                current_dict = Dict{String, Any}()
                in_boundary = false
            end
        end
    end
    
    # Verify we found all boundaries
    if n_boundaries > 0 && length(boundaries) != n_boundaries
        @warn "Expected $n_boundaries boundaries, found $(length(boundaries))"
    end
    
    return boundaries
end


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


