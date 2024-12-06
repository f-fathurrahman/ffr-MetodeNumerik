include("setup_works.jl")
import LinearAlgebra

include("material_point_type.jl")

function test_main()

    mpm_point = mpmMaterialPoint_2D_Classic()

    fOffset = 60.0/50/2.0
    thisMaterialDomain_01 = createMaterialDomain_Circle([30.0; 50.0], 9.6/2.0, fOffset);

    thisMaterialDomain_02 = createMaterialDomain_Rectangle([30.0; 20.0], 60.0, 40.6, fOffset);

end
test_main()