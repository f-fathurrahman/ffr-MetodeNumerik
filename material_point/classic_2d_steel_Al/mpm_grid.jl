
#fTime = 0.0 # global var?

function index2DTo1D(i::Int64, j::Int64, nColumns::Int64, nRows::Int64)
    index = nColumns*(j-1) + i
    if(index > nRows*nColumns || index < 1)
        @error("Index out of bounds: i=$(i), j=$(j)")
    end
    return index
end

mutable struct mpmGridPoint
    v2Fixed::Vector{Bool}
    fMass::Float64
    v2Position::Vector{Float64}
    v2Velocity::Vector{Float64}
    v2Momentum::Vector{Float64}
    v2Force::Vector{Float64}
end

function mpmGridPoint()
    v2Fixed = [false; false]
    fMass = 0.0
    v2Position = zeros(Float64, 2)
    v2Velocity = zeros(Float64, 2)
    v2Momentum = zeros(Float64, 2)
    v2Force = zeros(Float64, 2)
    return mpmGridPoint(
        v2Fixed, fMass, v2Position, v2Velocity, v2Momentum, v2Force
    )
end

mutable struct mpmGrid
    v2Length_Grid::Array{Float64}
    v2Nodes::Array{Int64}
    iNodes::Int64
    v2Length_Cell::Array{Float64}
    GridPoints::Vector{mpmGridPoint}
end

function mpmGrid(fGL_x::Float64, fGL_y::Float64, iN_x::Int64, iN_y::Int64)
    v2CL = zeros(Float64, 2)
    v2CL[1] = fGL_x / Float64(iN_x - 1.0)
    v2CL[2] = fGL_y / Float64(iN_y - 1.0)
    #
    thisGridPoint = Vector{mpmGridPoint}(undef, iN_x * iN_y)

    for j in 1:iN_y, i in 1:iN_x
        
        x = (i-1) * v2CL[1]
        y = (j-1) * v2CL[2]    
        index = index2DTo1D(i, j, iN_x, iN_y)

        bFixed_x = false
        bFixed_y = false
        
        if (i == 1) || (i == iN_x)
            bFixed_x = true
            bFixed_y = true
        end
        
        if (j == 1) || (j == iN_y)
            bFixed_x = true
            bFixed_y = true
        end

        thisGridPoint[index] = mpmGridPoint()
        thisGridPoint[index].v2Fixed[:] .= [bFixed_x; bFixed_y]
        thisGridPoint[index].v2Position[:] .= [x; y]
    end

    return mpmGrid(
        [fGL_x; fGL_y],
        [iN_x; iN_y],
        iN_x*iN_y,
        v2CL,
        thisGridPoint
    )
end

function getAdjacentGridPoints(
    thisMaterialPoint::mpmMaterialPoint_2D_Classic,
    thisGrid::mpmGrid
)
    v2Coordinate = thisMaterialPoint.v2Centroid
    fLength_Cell_x = thisGrid.v2Length_Cell[1]
    fLength_Cell_y = thisGrid.v2Length_Cell[2]

    iBottomLeft_i    = (floor(v2Coordinate[1] / fLength_Cell_x) + 1.0)
    iBottomLeft_j    = (floor(v2Coordinate[2] / fLength_Cell_y) + 1.0)

    if (iBottomLeft_j < 1) || (iBottomLeft_j > thisGrid.v2Nodes[2])
        error("Index out of bounds: j: $(iBottomLeft_j) v2Coordinate[2]: $(v2Coordinate[2])")
    end

    range1 = iBottomLeft_i:iBottomLeft_i+1
    range2 = iBottomLeft_j:iBottomLeft_j+1
    Nadjacent = length(range1) * length(range2)
    thisAdjacentGridPoints = Array{Int64}(undef, Nadjacent)
    ip = 0
    for i in range1, j in range2
        iIndex = index2DTo1D(Int64(i), Int64(j), thisGrid.v2Nodes[1], thisGrid.v2Nodes[2])
        ip += 1
        thisAdjacentGridPoints[ip] = iIndex
    end

    return thisAdjacentGridPoints
end
