@doc raw"""
    function quanticscrossinterpolate(
        ::Type{ValueType},
        f,
        grid::QuanticsGrids.Grid{n},
        initialpivots::Union{Nothing,AbstractVector{<:AbstractVector}}=nothing;
        nrandominitpivot=5,
        kwargs...
    ) where {ValueType}

Interpolate a function ``f(\mathbf{x})`` as a quantics tensor train. The tensor train itself is constructed using the 2-site tensor cross interpolation algorithm implemented in [`TensorCrossInterpolation.crossinterpolate2`](https://tensor4all.github.io/TensorCrossInterpolation.jl/dev/documentation/#TensorCrossInterpolation.crossinterpolate2-Union{Tuple{N},%20Tuple{ValueType},%20Tuple{Type{ValueType},%20Any,%20Union{Tuple{Vararg{Int64,%20N}},%20Vector{Int64}}},%20Tuple{Type{ValueType},%20Any,%20Union{Tuple{Vararg{Int64,%20N}},%20Vector{Int64}},%20Vector{Vector{Int64}}}}%20where%20{ValueType,%20N}).

Arguments:
- `ValueType` is the return type of `f`. Automatic inference is too error-prone.
- `f` is the function to be interpolated. `f` may take multiple arguments. The return type should be `ValueType`.
- `grid` is a `Grid{n}` object from [`QuanticsGrids.jl`](https://github.com/tensor4all/QuanticsGrids.jl) that describes a d-dimensional grid of discrete points indexed by binary digits. To avoid constructing a grid explicitly, use one of the other overloads.
- `initialpivots` is a vector of pivots to be used for initialization.
- `nrandominitpivot` determines how many random pivots should be used for initialization if no initial pivot is given.

All other arguments are forwareded to `crossinterpolate2`. Most importantly:
- `tolerance::Float64` is a float specifying the target tolerance for the interpolation. Default: `1e-8`.
- `pivottolerance::Float64` is a float that specifies the tolerance for adding new pivots, i.e. the truncation of tensor train bonds. It should be <= tolerance, otherwise convergence may be impossible. Default: `tolerance`.
- `maxbonddim::Int` specifies the maximum bond dimension for the TCI. Default: `typemax(Int)`, i.e. effectively unlimited.
- `maxiter::Int` is the maximum number of iterations (i.e. optimization sweeps) before aborting the TCI construction. Default: `200`.

For all other arguments, see the documentation for [`TensorCrossInterpolation.crossinterpolate2`](https://tensor4all.github.io/TensorCrossInterpolation.jl/dev/documentation/#TensorCrossInterpolation.crossinterpolate2-Union{Tuple{N},%20Tuple{ValueType},%20Tuple{Type{ValueType},%20Any,%20Union{Tuple{Vararg{Int64,%20N}},%20Vector{Int64}}},%20Tuple{Type{ValueType},%20Any,%20Union{Tuple{Vararg{Int64,%20N}},%20Vector{Int64}},%20Vector{Vector{Int64}}}}%20where%20{ValueType,%20N}).
"""



@doc raw"""
    function quanticscrossinterpolate(
        ::Type{ValueType},
        f,
        xvals::AbstractVector{<:AbstractVector},
        initialpivots::Union{Nothing,AbstractVector{<:AbstractVector}}=nothing;
        unfoldingscheme::Symbol=:interleaved,
        nrandominitpivot=5,
        kwargs...
    ) where {ValueType}

Interpolate a function ``f(\mathbf{x})`` as a quantics tensor train. This overload automatically constructs a Grid object from the ``\mathbf{x}`` points given in `xvals`.

Arguments:
- `xvals::AbstractVector{<:AbstractVector}`: A set of discrete points where `f` can be evaluated, given as a set of arrays, where `xvals[i]` describes the `i`th axis. Each array in `xvals` should contain `2^R` points for some integer `R`.
- For all other arguments, see the documentation of the main overload.
"""

