include("setup_works.jl")

using Printf
import LinearAlgebra

include("material_point_type.jl")
include("mpm_grid.jl")

thisGrid = mpmGrid(60.0, 60.0, 51, 51)

