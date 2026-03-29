using Pkg
Pkg.activate("QUANTICS", shared=true)

using Revise

const DIR_PACKAGES = ["./", "./MyTensorCrossInterpolation"]
for d in DIR_PACKAGES
    !( d in LOAD_PATH) && push!(LOAD_PATH, d)
end
