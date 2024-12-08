
function getShapeValue_Classic(
    thisMaterialPoint::mpmMaterialPoint_2D_Classic,
    thisGridPoint::MPMGridPoint,
    thisGrid::MPMGrid
)
    fShapeValue = 0.0

    v2Distance = thisMaterialPoint.centroid - thisGridPoint.position
    v2CellLength = thisGrid.cell_L;

    v2ShapeValue = zeros(2)
    v2ShapeValue[1] = 1.0 - abs(v2Distance[1]) / v2CellLength[1]
    v2ShapeValue[2] = 1.0 - abs(v2Distance[2]) / v2CellLength[2]

    if v2ShapeValue[1] < 0.0
        # @printf("Negative shape value!!! %e \n", v2Distance[1])
        v2ShapeValue[1] = 0.0
    end
    if v2ShapeValue[2] < 0.0
        # @printf("Negative shape value!!! \n")
        v2ShapeValue[2] = 0.0
    end

    fShapeValue = v2ShapeValue[1] * v2ShapeValue[2]
    #
    return fShapeValue
end

function getShapeGradient_Classic(
    thisMaterialPoint::mpmMaterialPoint_2D_Classic,
    thisGridPoint::MPMGridPoint,
    thisGrid::MPMGrid
)
    
    v2Result = zeros(2)

    v2Distance = thisMaterialPoint.centroid - thisGridPoint.position
    v2CellLength = thisGrid.cell_L;

    v2ShapeValue = zeros(2)
    v2ShapeValue[1] = 1.0 - abs(v2Distance[1]) / v2CellLength[1]
    v2ShapeValue[2] = 1.0 - abs(v2Distance[2]) / v2CellLength[2]

    if(v2ShapeValue[1] < 0.0)
        v2ShapeValue[1] = 0.0
    end
    if(v2ShapeValue[2] < 0.0)
        v2ShapeValue[2] = 0.0
    end

    v2Result[1] = -v2ShapeValue[2]*sign(v2Distance[1]) / v2CellLength[1]
    v2Result[2] = -v2ShapeValue[1]*sign(v2Distance[2]) / v2CellLength[2]

    return v2Result
end

