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