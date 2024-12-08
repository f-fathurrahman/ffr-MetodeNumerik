
#fTime = 0.0 # global var?

function index2DTo1D(i::Int64, j::Int64, nColumns::Int64, nRows::Int64)
    index = nColumns*(j-1) + i
    if(index > nRows*nColumns || index < 1)
        @error("Index out of bounds: i=$(i), j=$(j)")
    end
    return index
end

#
# FIXME: comibine MPMGridPoint with MPMGrid ?
#

mutable struct MPMGridPoint
    is_fixed::Vector{Bool}
    mass::Float64
    position::Vector{Float64}
    velocity::Vector{Float64}
    momentum::Vector{Float64}
    force::Vector{Float64}
end

function MPMGridPoint()
    is_fixed = [false; false]
    mass = 0.0
    position = zeros(Float64, 2)
    velocity = zeros(Float64, 2)
    momentum = zeros(Float64, 2)
    force = zeros(Float64, 2)
    return MPMGridPoint(
        is_fixed, mass, position, velocity, momentum, force
    )
end

mutable struct MPMGrid
    L::Array{Float64} # lengths
    Nnodes::Array{Int64}
    NnodesTotal::Int64
    cell_L::Array{Float64}
    points::Vector{MPMGridPoint}
end

function MPMGrid(
    grid_Lx::Float64,
    grid_Ly::Float64,
    Nnodes_x::Int64,
    Nnodes_y::Int64
)
    cell_L = zeros(Float64, 2)
    cell_L[1] = grid_Lx / Float64(Nnodes_x - 1.0)
    cell_L[2] = grid_Ly / Float64(Nnodes_y - 1.0)
    #
    points = Vector{MPMGridPoint}(undef, Nnodes_x * Nnodes_y)
    for j in 1:Nnodes_y, i in 1:Nnodes_x
        x = (i-1) * cell_L[1]
        y = (j-1) * cell_L[2]
        ip = index2DTo1D(i, j, Nnodes_x, Nnodes_y) # need this?
        #
        points[ip] = MPMGridPoint()
        points[ip].position[1] = x
        points[ip].position[2] = y
        #
        if (i == 1) || (i == Nnodes_x)
            points[ip].is_fixed[1] = true
            points[ip].is_fixed[2] = true
        end
        #
        if (j == 1) || (j == Nnodes_y)
            points[ip].is_fixed[1] = true
            points[ip].is_fixed[2] = true
        end
    end

    return MPMGrid(
        [grid_Lx; grid_Ly],
        [Nnodes_x; Nnodes_y],
        Nnodes_x*Nnodes_y,
        cell_L,
        points
    )
end

function getAdjacentGridPoints(
    thisMaterialPoint::mpmMaterialPoint_2D_Classic,
    thisGrid::MPMGrid
)
    v2Coordinate = thisMaterialPoint.centroid
    fLength_Cell_x = thisGrid.cell_L[1]
    fLength_Cell_y = thisGrid.cell_L[2]

    iBottomLeft_i    = (floor(v2Coordinate[1] / fLength_Cell_x) + 1.0)
    iBottomLeft_j    = (floor(v2Coordinate[2] / fLength_Cell_y) + 1.0)

    if (iBottomLeft_j < 1) || (iBottomLeft_j > thisGrid.Nnodes[2])
        error("Index out of bounds: j: $(iBottomLeft_j) v2Coordinate[2]: $(v2Coordinate[2])")
    end

    range1 = iBottomLeft_i:iBottomLeft_i+1
    range2 = iBottomLeft_j:iBottomLeft_j+1
    Nadjacent = length(range1) * length(range2)
    thisAdjacentGridPoints = Array{Int64}(undef, Nadjacent)
    ip = 0
    for i in range1, j in range2
        iIndex = index2DTo1D(Int64(i), Int64(j), thisGrid.Nnodes[1], thisGrid.Nnodes[2])
        ip += 1
        thisAdjacentGridPoints[ip] = iIndex
    end

    return thisAdjacentGridPoints
end
