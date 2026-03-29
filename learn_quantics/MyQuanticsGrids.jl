module MyQuanticsGrids

export DiscretizedGrid, InherentDiscreteGrid
export quantics_to_grididx, quantics_to_origcoord
export grididx_to_quantics, grididx_to_origcoord
export origcoord_to_quantics, origcoord_to_grididx

abstract type Grid{D} end

struct InherentDiscreteGrid{D} <: Grid{D}
    Rs::NTuple{D,Int}
    origin::NTuple{D,Int}
    step::NTuple{D,Int}
    variablenames::NTuple{D,Symbol}
    base::Int
    indextable::Vector{Vector{Tuple{Symbol,Int}}}

    # Lookup table: lookup_table[variablename_index][bitnumber] -> (site_index, position_in_site)
    lookup_table::NTuple{D,Vector{Tuple{Int,Int}}}
    maxgrididx::NTuple{D,Int}

    function InherentDiscreteGrid{D}(
        Rs::NTuple{D,Int},
        origin::NTuple{D,Int},
        step::NTuple{D,Int},
        variablenames::NTuple{D,Symbol},
        base::Int,
        indextable::Vector{Vector{Tuple{Symbol,Int}}}
    ) where {D}
        if !(D isa Int)
            throw(ArgumentError(lazy"Got dimension $D, which is not an Int."))
        end
        if base <= 1
            throw(ArgumentError(lazy"Got base = $base. base must be at least 2."))
        end
        for d in 1:D
            if !(step[d] >= 1)
                throw(ArgumentError(lazy"Got step[$d] = $(step[d]). step must be at least 1."))
            end
            if !rangecheck_R(Rs[d]; base)
                throw(ArgumentError(lazy"Got Rs[$d] = $(Rs[d]) with base = $base. For all gridindices from 1 to base^R to fit into an Int, we must have base ^ R <= typemax(Int)"))
            end
        end
        if !allunique(variablenames)
            throw(ArgumentError(lazy"Got variablenames = $variablenames. variablenames must be unique."))
        end

        lookup_table = _build_lookup_table(Rs, indextable, variablenames)
        maxgrididx = map(R -> base^R, Rs)

        return new{D}(Rs, origin, step, variablenames, base, indextable, lookup_table, maxgrididx)
    end
end

# ============================================================================
# Helper/utility functions
# ============================================================================

_convert_to_scalar_if_possible(x) = x
_convert_to_scalar_if_possible(x::NTuple{1,T}) where {T} = first(x)

_to_tuple(::Val{d}, x::NTuple{d}) where {d} = x
_to_tuple(::Val{d}, x) where {d} = ntuple(i -> x, d)

default_step(::Val{D}) where D = ntuple(Returns(1), D)

default_origin(::Val{D}) where D = ntuple(Returns(1), D)

function _build_lookup_table(Rs::NTuple{D,Int}, indextable::Vector{Vector{Tuple{Symbol,Int}}}, variablenames::NTuple{D,Symbol}) where D
    lookup_table = ntuple(D) do d
        if Rs[d] < 0
            throw(ArgumentError(lazy"Got Rs[$d] = $(Rs[d]). Rs must be non-negative."))
        end
        Vector{Tuple{Int,Int}}(undef, Rs[d])
    end

    index_visited = [fill(false, Rs[d]) for d in 1:D]

    for (site_idx, quanticsindices) in pairs(indextable)
        for (pos_in_site, qindex) in pairs(quanticsindices)
            variablename, bitnumber = qindex
            var_idx = findfirst(==(variablename), variablenames)
            if isnothing(var_idx)
                throw(ArgumentError(lazy"Index table contains unknown index $qindex. Valid variablenames are $variablenames."))
            elseif bitnumber > Rs[var_idx]
                throw(ArgumentError(lazy"Index table contains quantics index $bitnumber of variable $variablename, but it must be smaller than or equal to the number of quantics indices for that variable, which is $(Rs[var_idx])."))
            elseif index_visited[var_idx][bitnumber]
                throw(ArgumentError(lazy"Index table contains quantics index $bitnumber of variable $variablename more than once."))
            end

            lookup_table[var_idx][bitnumber] = (site_idx, pos_in_site)
            index_visited[var_idx][bitnumber] = true
        end
    end

    for (var_idx, visited) in enumerate(index_visited)
        bitnumber = findfirst(==(false), visited)
        if !isnothing(bitnumber)
            throw(ArgumentError(lazy"Index table contains no site for quantics index $bitnumber of variable $(variablenames[var_idx])."))
        end
    end

    return lookup_table
end

function rangecheck_R(R::Int; base::Int=2)::Bool
    # For all gridindices from 1 to base^R to fit into an Int64,
    # we must have base ^ R <= typemax(Int)
    result = 1
    for _ in 1:R
        result <= typemax(Int) ÷ base || return false
        result *= base
    end
    return true
end

function _build_indextable(variablenames::NTuple{D,Symbol}, Rs::NTuple{D,Int}, unfoldingscheme::Symbol) where D
    if !(unfoldingscheme in (:interleaved, :fused))
        throw(ArgumentError(lazy"Got unfoldingscheme = $unfoldingscheme. Supported are :interleaved and :fused."))
    end
    indextable = Vector{Tuple{Symbol,Int}}[]

    for bitnumber in 1:maximum(Rs; init=0)
        if unfoldingscheme === :interleaved
            _add_interleaved_indices!(indextable, variablenames, Rs, bitnumber)
        elseif unfoldingscheme === :fused
            _add_fused_indices!(indextable, variablenames, Rs, bitnumber)
        end
    end

    return indextable
end

function _add_interleaved_indices!(indextable, variablenames::NTuple{D,Symbol}, Rs::NTuple{D,Int}, bitnumber) where D
    for d in 1:D
        bitnumber ∈ 1:Rs[d] || continue
        qindex = (variablenames[d], bitnumber)
        push!(indextable, [qindex])
    end
end

function _add_fused_indices!(indextable, variablenames::NTuple{D,Symbol}, Rs::NTuple{D,Int}, bitnumber) where D
    indices_bitnumber = Tuple{Symbol,Int}[]
    # Add dimensions in reverse order to match DiscretizedGrid convention
    # where the first dimension varies fastest in fused quantics
    for d in D:-1:1
        bitnumber ∈ 1:Rs[d] || continue
        qindex = (variablenames[d], bitnumber)
        push!(indices_bitnumber, qindex)
    end
    if !isempty(indices_bitnumber)
        push!(indextable, indices_bitnumber)
    end
end

function _quantics_to_grididx_general(g::InherentDiscreteGrid{D}, quantics) where D
    base = g.base

    return ntuple(D) do d
        grididx = 1
        R_d = g.Rs[d]

        for bitnumber in 1:R_d
            site_idx, pos_in_site = g.lookup_table[d][bitnumber]
            quantics_val = quantics[site_idx]
            site_len = length(g.indextable[site_idx])

            temp = quantics_val - 1
            for _ in 1:(site_len-pos_in_site)
                temp = div(temp, base)
            end
            digit = temp % base

            grididx += digit * base^(R_d - bitnumber)
        end
        grididx
    end
end

function _quantics_to_grididx_base2(g::InherentDiscreteGrid{D}, quantics) where D
    return ntuple(D) do d
        grididx = 0
        R_d = g.Rs[d]

        for bitnumber in 1:R_d
            site_idx, pos_in_site = g.lookup_table[d][bitnumber]
            bit_position = length(g.indextable[site_idx]) - pos_in_site
            digit = ((quantics[site_idx] - 1) >> bit_position) & 1
            grididx |= digit << (R_d - bitnumber)
        end
        grididx + 1
    end
end

function _grididx_to_quantics_general!(result::Vector{Int}, g::InherentDiscreteGrid{D}, grididx::NTuple{D,Int}) where D
    base = g.base

    @inbounds for d in 1:D
        zero_based_idx = grididx[d] - 1
        R_d = g.Rs[d]

        for bitnumber in 1:R_d
            site_idx, pos_in_site = g.lookup_table[d][bitnumber]
            site_length = length(g.indextable[site_idx])

            bit_position = R_d - bitnumber
            digit = (zero_based_idx ÷ (base^bit_position)) % base

            power = site_length - pos_in_site
            result[site_idx] += digit * (base^power)
        end
    end
end

function _grididx_to_quantics_base2!(result::Vector{Int}, g::InherentDiscreteGrid{D}, grididx::NTuple{D,Int}) where D
    @inbounds for d in 1:D
        zero_based_idx = grididx[d] - 1
        R_d = g.Rs[d]

        for bitnumber in 1:R_d
            site_idx, pos_in_site = g.lookup_table[d][bitnumber]
            site_length = length(g.indextable[site_idx])

            bit_position = R_d - bitnumber
            digit = (zero_based_idx >> bit_position) & 1

            power = site_length - pos_in_site
            result[site_idx] += digit << power
        end
    end
end

# ============================================================================
# Constructors
# ============================================================================

function InherentDiscreteGrid{D}(
    Rs,
    origin=default_origin(Val(D));
    unfoldingscheme=:fused,
    step=default_step(Val(D)),
    base=2
) where D
    Rs = _to_tuple(Val(D), Rs)
    origin = _to_tuple(Val(D), origin)
    step = _to_tuple(Val(D), step)
    variablenames = ntuple(Symbol, D)
    indextable = _build_indextable(variablenames, Rs, unfoldingscheme)
    InherentDiscreteGrid{D}(Rs, origin, step, variablenames, base, indextable)
end

function InherentDiscreteGrid(
    variablenames::NTuple{D,Symbol},
    indextable::Vector{Vector{Tuple{Symbol,Int}}};
    origin=default_origin(Val(D)),
    step=default_step(Val(D)),
    base=2
) where D
    Rs = Tuple(map(variablenames) do variablename
        count(index -> first(index) == variablename, Iterators.flatten(indextable))
    end)

    return InherentDiscreteGrid{D}(Rs, origin, step, variablenames, base, indextable)
end

function InherentDiscreteGrid(
    variablenames::NTuple{D,Symbol},
    Rs::NTuple{D,Int};
    origin=default_origin(Val(D)),
    step=default_step(Val(D)),
    base=2,
    unfoldingscheme=:fused
) where {D}
    origin = _to_tuple(Val(D), origin)
    step = _to_tuple(Val(D), step)
    indextable = _build_indextable(variablenames, Rs, unfoldingscheme)
    return InherentDiscreteGrid{D}(Rs, origin, step, variablenames, base, indextable)
end

function InherentDiscreteGrid(Rs::NTuple{D,Int}; variablenames=ntuple(Symbol, D), kwargs...) where {D}
    return InherentDiscreteGrid(variablenames, Rs; kwargs...)
end

function InherentDiscreteGrid(R::Int, origin; kwargs...)
    return InherentDiscreteGrid{length(origin)}(R, origin; kwargs...)
end

# ============================================================================
# Basic property accessor functions
# ============================================================================

Base.ndims(::InherentDiscreteGrid{D}) where D = D

Base.length(g::InherentDiscreteGrid) = length(g.indextable)

grid_Rs(g::InherentDiscreteGrid) = g.Rs

grid_indextable(g::InherentDiscreteGrid) = g.indextable

grid_base(g::InherentDiscreteGrid) = g.base

grid_variablenames(g::InherentDiscreteGrid) = g.variablenames

grid_step(g::InherentDiscreteGrid) = _convert_to_scalar_if_possible(g.step)

grid_origin(g::InherentDiscreteGrid) = _convert_to_scalar_if_possible(g.origin)

function sitedim(g::InherentDiscreteGrid, site::Int)::Int
    if !(site ∈ eachindex(g.indextable))
        throw(DomainError(site, lazy"Site index out of bounds [1, $(length(g.indextable))]."))
    end
    return g.base^length(g.indextable[site])
end

# ============================================================================
# Grid coordinate functions
# ============================================================================

grid_min(g::InherentDiscreteGrid) = _convert_to_scalar_if_possible(g.origin)

grid_max(g::InherentDiscreteGrid) = _convert_to_scalar_if_possible(
    grid_origin(g) .+ grid_step(g) .* (g.base .^ g.Rs .- 1),
)

# ============================================================================
# Core conversion functions
# ============================================================================

function quantics_to_grididx(g::InherentDiscreteGrid{D}, quantics::AbstractVector{Int}) where D
    # TODO: add switch to turn off input validation
    if !(length(quantics) == length(g))
        throw(ArgumentError(lazy"Quantics vector must have length $(length(g.indextable)), got $(length(quantics))."))
    end

    for site in eachindex(quantics)
        if !(1 <= quantics[site] <= sitedim(g, site))
            throw(DomainError(quantics[site], lazy"Quantics value for site $site out of range 1:$(sitedim(g, site))."))
        end
    end

    result = if g.base == 2
        _quantics_to_grididx_base2(g, quantics)
    else
        _quantics_to_grididx_general(g, quantics)
    end

    return _convert_to_scalar_if_possible(result)
end

function grididx_to_quantics(g::InherentDiscreteGrid{D}, grididx::Int) where D
    grididx_tuple = _to_tuple(Val(D), grididx)
    return grididx_to_quantics(g, grididx_tuple)
end
function grididx_to_quantics(g::InherentDiscreteGrid{D}, grididx_tuple::NTuple{D,Int}) where D
    # TODO: add switch to turn off input validation
    for d in 1:D
        if !(1 <= grididx_tuple[d] <= g.maxgrididx[d])
            throw(DomainError(grididx_tuple[d], lazy"Grid index out of bounds [1, $(g.maxgrididx[d])]."))
        end
    end

    result = ones(Int, length(g.indextable))
    if g.base == 2
        _grididx_to_quantics_base2!(result, g, grididx_tuple)
    else
        _grididx_to_quantics_general!(result, g, grididx_tuple)
    end
    return result
end

function grididx_to_origcoord(g::InherentDiscreteGrid{D}, grididx) where D
    grididx = _to_tuple(Val(D), grididx)
    for d in 1:D
        if !(grididx[d] ∈ 1:g.maxgrididx[d])
            throw(DomainError(grididx[d], lazy"Grid index out of bounds [1, $(g.maxgrididx[d])]."))
        end
    end

    res = grid_origin(g) .+ (grididx .- 1) .* grid_step(g)

    return _convert_to_scalar_if_possible(res)
end

function origcoord_to_grididx(g::InherentDiscreteGrid{D}, coordinate) where {D}
    coord_tuple = _to_tuple(Val(D), coordinate)
    bounds_lower = grid_min(g)
    bounds_upper = grid_max(g)
    # TODO: think about the correct bounds to use here
    for d in 1:D
        if !(bounds_lower[d] <= coord_tuple[d] <= bounds_upper[d])
            throw(DomainError(coord_tuple[d], lazy"Coordinate out of bounds [$(bounds_lower[d]), $(bounds_upper[d])]."))
        end
    end

    discrete_idx = div.(coord_tuple .- grid_min(g), grid_step(g)) .+ 1
    discrete_idx = clamp.(discrete_idx, 1, g.base .^ g.Rs)

    return _convert_to_scalar_if_possible(discrete_idx)
end

function origcoord_to_quantics(g::InherentDiscreteGrid, coordinate)
    grididx_to_quantics(g, origcoord_to_grididx(g, coordinate))
end

function quantics_to_origcoord(g::InherentDiscreteGrid, quantics)
    grididx_to_origcoord(g, quantics_to_grididx(g, quantics))
end

# ============================================================================
# Other utility functions
# ============================================================================

function localdimensions(g::InherentDiscreteGrid)::Vector{Int}
    return g.base .^ length.(g.indextable)
end


"""
    DiscretizedGrid{D}

A discretized grid structure for D-dimensional grids with variable resolution,
supporting efficient conversion between quantics, grid indices, and original coordinates.
A `DiscretizedGrid` instance is intended to undergird a quantics tensor train
with a specific index structure, as defined in the `indextable` field.
For example, say indextable is `[[(:a, 1), (:b, 2)], [(:a, 2)], [(:b, 1), (:a, 3)]]`,
then the corresponding tensor train has 3 tensor cores:

      a_1 b_2          a_2          b_1 a_3
       |   |            |            |   |
    ┌──┴───┴──┐    ┌────┴────┐    ┌──┴───┴──┐
    │         │    │         │    │         │
    │         │────│         │────│         │
    │         │    │         │    │         │
    └─────────┘    └─────────┘    └─────────┘

This object may be constructed with
```julia-repl
julia> grid = DiscretizedGrid((:a, :b), [[(:a, 1), (:b, 2)], [(:a, 2)], [(:b, 1), (:a, 3)]])
DiscretizedGrid{2} with 8×4 = 32 grid points
├─ Variables: (a, b)
├─ Resolutions: (a: 3, b: 2)
├─ Domain: unit square [0, 1)²
├─ Grid spacing: (Δa = 0.125, Δb = 0.25)
└─ Tensor train: 3 sites (dimensions: 4-2-4)
```
and represents a 2^3 x 2^2 discretization of the unit square in the 2D plane (the x mark grid points):

       1.0  ┌───────────────────────────────┐
            │                               │
            │                               │
       0.75 x   x   x   x   x   x   x   x   │
            │                               │
            │                               │
    b  0.5  x   x   x   x   x   x   x   x   │
            │                               │
            │                               │
       0.25 x   x   x   x   x   x   x   x   │
            │                               │
            │                               │
       0.0  x───x───x───x───x───x───x───x───┘
           0.0     0.25    0.5     0.75    1.0

                            a

If something other than a unit square is desired, `lower_bound` and `upper_bound`
can be specified. Also, bases different than the default base 2 can be used,
e.g. `base=3` for a ternary grid.

In addition to the plain methods, there is a convenience layer for conversion
from the original coordinates
```julia-repl
julia> origcoord_to_grididx(grid; a=0.5, b=0.25)
(5, 2)

julia> origcoord_to_quantics(grid; a=0.5, b=0.25)
3-element Vector{Int64}:
 4
 1
 1
```
and also from grid indices
```julia-repl
julia> grididx_to_origcoord(grid; a=5, b=2)
(0.5, 0.25)

julia> grididx_to_quantics(grid; a=5, b=2)
3-element Vector{Int64}:
 4
 1
 1
```

For a simpler grid, we can just supply the resolution in each dimension:
```julia-repl
julia> boring_grid = DiscretizedGrid((3, 9))
DiscretizedGrid{2} with 8×512 = 4096 grid points
├─ Resolutions: (1: 3, 2: 9)
├─ Domain: unit square [0, 1)²
├─ Grid spacing: (Δ1 = 0.125, Δ2 = 0.001953125)
└─ Tensor train: 12 sites (uniform dimension 2)
```
In this case, variable names are automatically generated as `1`, `2`, etc.
"""
struct DiscretizedGrid{D} <: Grid{D}
    discretegrid::InherentDiscreteGrid{D}
    lower_bound::NTuple{D,Float64}
    upper_bound::NTuple{D,Float64}

    function DiscretizedGrid{D}(
        Rs, lower_bound, upper_bound, variablenames, base, indextable, includeendpoint
    ) where {D}
        lower_bound = _to_tuple(Val(D), lower_bound)
        upper_bound = _to_tuple(Val(D), upper_bound)
        for d in 1:D
            if !(lower_bound[d] < upper_bound[d])
                throw(ArgumentError(lazy"Got (lower_bound[$d], upper_bound[$d]) = $((lower_bound[d], upper_bound[d])). Each lower bound needs to be strictly less than the corresponding upper bound."))
            end
        end

        discretegrid = InherentDiscreteGrid(variablenames, indextable; base)

        upper_bound = _adjust_upper_bounds(
            upper_bound, lower_bound, includeendpoint, base, Rs, Val(D)
        )

        return new{D}(discretegrid, lower_bound, upper_bound)
    end
end

# ============================================================================
# Helper/utility functions
# ============================================================================

function _adjust_upper_bounds(upper_bound, lower_bound, includeendpoint, base, Rs, ::Val{D}) where D
    includeendpoint = _to_tuple(Val(D), includeendpoint)
    for d in 1:D
        if iszero(Rs[d]) && includeendpoint[d]
            throw(ArgumentError(lazy"Got Rs[$d] = 0 and includeendpoint[$d] = true. This is not allowed."))
        end
    end

    return ntuple(D) do d
        if includeendpoint[d]
            upper_bound[d] + (upper_bound[d] - lower_bound[d]) / (base^Rs[d] - 1)
        else
            upper_bound[d]
        end
    end
end

default_lower_bound(::Val{D}) where D = _to_tuple(Val(D), 0.0)

default_upper_bound(::Val{D}) where D = _to_tuple(Val(D), 1.0)

function _handle_kwargs_input(g::DiscretizedGrid{D}; kwargs...) where {D}
    provided_keys = keys(kwargs)
    expected_keys = grid_variablenames(g)
    if !(Set(provided_keys) == Set(expected_keys))
        throw(ArgumentError(lazy"Expected keyword arguments $(expected_keys), got $(tuple(provided_keys...))"))
    end

    if !all(v -> v isa Real, values(kwargs))
        throw(ArgumentError(lazy"Got kwargs = $kwargs. All keyword argument values must be Real numbers"))
    end

    return ntuple(D) do d
        variablename = grid_variablenames(g)[d]
        kwargs[variablename]
    end
end

# ============================================================================
# Constructors
# ============================================================================

function DiscretizedGrid(
    variablenames::NTuple{D,Symbol},
    Rs::NTuple{D,Int};
    lower_bound=default_lower_bound(Val(D)),
    upper_bound=default_upper_bound(Val(D)),
    base::Int=2,
    unfoldingscheme::Symbol=:fused,
    includeendpoint=false
) where {D}
    indextable = _build_indextable(variablenames, Rs, unfoldingscheme)

    return DiscretizedGrid{D}(Rs, lower_bound, upper_bound, variablenames, base, indextable, includeendpoint)
end

function DiscretizedGrid(Rs::NTuple{D,Int}; variablenames=ntuple(Symbol, D), kwargs...) where {D}
    return DiscretizedGrid(variablenames, Rs; kwargs...)
end

function DiscretizedGrid{D}(
    R,
    lower_bound=default_lower_bound(Val(D)),
    upper_bound=default_upper_bound(Val(D));
    kwargs...
) where {D}
    return DiscretizedGrid(_to_tuple(Val(D), R); lower_bound, upper_bound, kwargs...)
end

function DiscretizedGrid(
    R,
    lower_bound::NTuple{D,Real},
    upper_bound::NTuple{D,Real};
    kwargs...
) where {D}
    return DiscretizedGrid(_to_tuple(Val(D), R); lower_bound, upper_bound, kwargs...)
end

function DiscretizedGrid(
    variablenames::NTuple{D,Symbol},
    indextable::Vector{Vector{Tuple{Symbol,Int}}};
    lower_bound=default_lower_bound(Val(D)),
    upper_bound=default_upper_bound(Val(D)),
    base::Int=2,
    includeendpoint=false
) where D
    Rs = Tuple(map(variablenames) do variablename
        count(index -> first(index) == variablename, Iterators.flatten(indextable))
    end)

    return DiscretizedGrid{D}(Rs, lower_bound, upper_bound, variablenames, base, indextable, includeendpoint)
end

function DiscretizedGrid(R::Int, lower_bound::Real, upper_bound::Real; kwargs...)
    return DiscretizedGrid{1}(R, lower_bound, upper_bound; kwargs...)
end

# ============================================================================
# Basic property accessor functions
# ============================================================================

Base.ndims(::DiscretizedGrid{D}) where D = D

Base.length(g::DiscretizedGrid) = length(grid_indextable(g))

grid_Rs(g::DiscretizedGrid{D}) where D = grid_Rs(g.discretegrid)

grid_indextable(g::DiscretizedGrid{D}) where D = grid_indextable(g.discretegrid)

grid_base(g::DiscretizedGrid{D}) where D = grid_base(g.discretegrid)

grid_variablenames(g::DiscretizedGrid{D}) where D = grid_variablenames(g.discretegrid)

upper_bound(g::DiscretizedGrid) = _convert_to_scalar_if_possible(g.upper_bound)

lower_bound(g::DiscretizedGrid) = _convert_to_scalar_if_possible(g.lower_bound)

grid_origin(g::DiscretizedGrid) = lower_bound(g)

function sitedim(g::DiscretizedGrid, site::Int)::Int
    if !(site ∈ eachindex(grid_indextable(g)))
        throw(DomainError(site, lazy"Site index out of bounds [1, $(length(grid_indextable(g)))]"))
    end
    return grid_base(g)^length(grid_indextable(g)[site])
end

# ============================================================================
# Grid coordinate functions
# ============================================================================

grid_min(g::DiscretizedGrid) = _convert_to_scalar_if_possible(g.lower_bound)

grid_max(g::DiscretizedGrid) = _convert_to_scalar_if_possible(g.upper_bound .- grid_step(g))

grid_step(g::DiscretizedGrid) = _convert_to_scalar_if_possible(
    (upper_bound(g) .- lower_bound(g)) ./ (grid_base(g) .^ grid_Rs(g)),
)

function grid_origcoords(g::DiscretizedGrid, d::Int)
    if !(1 ≤ d ≤ ndims(g))
        throw(DomainError(d, lazy"Dimension $d out of bounds [1, $(ndims(g))]."))
    end
    start = grid_min(g)[d]
    stop = grid_max(g)[d]
    length = grid_base(g)^grid_Rs(g)[d]
    return range(start, stop, length)
end

function grid_origcoords(g::DiscretizedGrid, variablename::Symbol)
    d = findfirst(==(variablename), grid_variablenames(g))
    isnothing(d) && throw(ArgumentError(lazy"Variable name :$variablename not found in grid. Available variables: $(grid_variablenames(g))"))
    return grid_origcoords(g, d)
end

# ============================================================================
# Core conversion functions
# ============================================================================

function quantics_to_grididx(g::DiscretizedGrid, quantics::AbstractVector{Int})
    return quantics_to_grididx(g.discretegrid, quantics)
end

function grididx_to_quantics(g::DiscretizedGrid, grididx)
    return grididx_to_quantics(g.discretegrid, grididx)
end

function grididx_to_origcoord(g::DiscretizedGrid{D}, index) where {D}
    index_tuple = _to_tuple(Val(D), index)
    for d in 1:D
        if !(1 <= index_tuple[d] <= grid_base(g)^grid_Rs(g)[d])
            throw(DomainError(index_tuple[d], lazy"Grid index out of bounds [1, $(grid_base(g) ^ grid_Rs(g)[d])]."))
        end
    end

    res = ntuple(D) do d
        step_d = (g.upper_bound[d] - g.lower_bound[d]) / (grid_base(g)^grid_Rs(g)[d])
        g.lower_bound[d] + (index_tuple[d] - 1) * step_d
    end

    return _convert_to_scalar_if_possible(res)
end

function origcoord_to_grididx(g::DiscretizedGrid{D}, coordinate) where {D}
    coord_tuple = _to_tuple(Val(D), coordinate)

    bounds_lower = lower_bound(g)
    bounds_upper = upper_bound(g)
    for d in 1:D
        if !(bounds_lower[d] <= coord_tuple[d] <= bounds_upper[d])
            throw(DomainError(coord_tuple[d], lazy"Coordinate out of bounds [$(bounds_lower[d]), $(bounds_upper[d])] (dimension $d)."))
        end
    end

    steps = grid_step(g)
    indices = ntuple(D) do d
        target = coord_tuple[d]
        step_d = steps[d]

        continuous_idx = (target - bounds_lower[d]) / step_d + 1

        discrete_idx = round(Int, continuous_idx)
        clamp(discrete_idx, 1, grid_base(g)^grid_Rs(g)[d])
    end

    return _convert_to_scalar_if_possible(indices)
end

function origcoord_to_quantics(g::DiscretizedGrid{D}, coordinate) where {D}
    coord_tuple = _to_tuple(Val(D), coordinate)
    grid_idx = origcoord_to_grididx(g, coord_tuple)
    return grididx_to_quantics(g, grid_idx)
end

function quantics_to_origcoord(g::DiscretizedGrid{D}, quantics::AbstractVector{Int}) where {D}
    grid_idx = quantics_to_grididx(g, quantics)
    return grididx_to_origcoord(g, grid_idx)
end

# ============================================================================
# Keyword argument convenience functions
# ============================================================================

function origcoord_to_grididx(g::DiscretizedGrid; kwargs...)
    coordinate = _handle_kwargs_input(g; kwargs...)
    return origcoord_to_grididx(g, coordinate)
end

function origcoord_to_quantics(g::DiscretizedGrid; kwargs...)
    coordinate = _handle_kwargs_input(g; kwargs...)
    return origcoord_to_quantics(g, coordinate)
end

function grididx_to_origcoord(g::DiscretizedGrid; kwargs...)
    index = _handle_kwargs_input(g; kwargs...)
    return grididx_to_origcoord(g, index)
end

function grididx_to_quantics(g::DiscretizedGrid; kwargs...)
    index = _handle_kwargs_input(g; kwargs...)
    return grididx_to_quantics(g, index)
end

# ============================================================================
# Other utility functions
# ============================================================================

function localdimensions(g::DiscretizedGrid)::Vector{Int}
    return grid_base(g) .^ length.(grid_indextable(g))
end

function quanticsfunction(::Type{T}, g::DiscretizedGrid, f::F)::Function where {T,F<:Function}
    function wrapped_function(quantics)::T
        coords = quantics_to_origcoord(g, quantics)
        if coords isa Tuple
            return f(coords...)
        else
            return f(coords)
        end
    end
    return wrapped_function
end

# ============================================================================
# Display/show methods
# ============================================================================

function Base.show(io::IO, ::MIME"text/plain", g::DiscretizedGrid{D}) where D
    print(io, "DiscretizedGrid{$D}")

    # Grid resolution and total points
    total_points = prod(big(grid_base(g)) .^ grid_Rs(g))
    if D <= 1
        print(io, " with $total_points grid point" * (isone(total_points) ? "" : "s"))
    else
        print(io, " with $(join(big(grid_base(g)) .^ grid_Rs(g), " × ")) = $total_points grid points")
    end

    # Variable names (if meaningful)
    if any(name -> !startswith(string(name), r"^\d+$"), grid_variablenames(g))
        var_str = join(grid_variablenames(g), ", ")
        print(io, "\n├─ Variables: ($var_str)")
    end

    # Resolution per dimension
    if D == 1
        print(io, "\n├─ Resolution: $(grid_Rs(g)[1]) bits")
    else
        res_str = join(["$(grid_variablenames(g)[i]): $(grid_Rs(g)[i])" for i in 1:D], ", ")
        print(io, "\n├─ Resolutions: ($res_str)")
    end

    # Bounds (only show if not default unit interval/square/cube)
    default_lower = default_lower_bound(Val(D))
    default_upper = default_upper_bound(Val(D))
    if lower_bound(g) != default_lower || any(abs.(upper_bound(g) .- default_upper) .> 1e-10)
        if D == 1
            print(io, "\n├─ Domain: [$(lower_bound(g)[1]), $(upper_bound(g)[1]))")
        else
            bounds_str = join(["[$(lower_bound(g)[i]), $(upper_bound(g)[i]))" for i in 1:D], " × ")
            print(io, "\n├─ Domain: $bounds_str")
        end

        # Grid spacing
        step_vals = grid_step(g)
        if D == 1
            print(io, "\n├─ Grid spacing: $(step_vals)")
        else
            step_str = join(["Δ$(grid_variablenames(g)[i]) = $(step_vals[i])" for i in 1:D], ", ")
            print(io, "\n├─ Grid spacing: ($step_str)")
        end
    else
        # For unit domain, show appropriate unit description
        unit_domain_str = if D == 1
            "unit interval [0, 1)"
        elseif D == 2
            "unit square [0, 1)²"
        elseif D == 3
            "unit cube [0, 1)³"
        else
            "unit hypercube [0, 1)^$D"
        end
        print(io, "\n├─ Domain: $unit_domain_str")

        # Grid spacing
        step_vals = grid_step(g)
        if D == 1
            print(io, "\n├─ Grid spacing: $(step_vals)")
        else
            step_str = join(["Δ$(grid_variablenames(g)[i]) = $(step_vals[i])" for i in 1:D], ", ")
            print(io, "\n├─ Grid spacing: ($step_str)")
        end
    end

    # Base (only show if not binary)
    if grid_base(g) != 2
        print(io, "\n├─ Base: $(grid_base(g))")
    end

    # Tensor structure summary
    num_sites = length(grid_indextable(g))
    sitedims = Int[sitedim(g, site) for site in 1:num_sites]

    print(io, "\n└─ Tensor train: $num_sites sites")
    if !isempty(sitedims)
        if allequal(sitedims)
            print(io, " (uniform dimension $(sitedims[1]))")
        else
            print(io, " (dimensions: $(join(sitedims, "-")))")
        end
    end
end

function Base.show(io::IO, g::DiscretizedGrid{D}) where D
    total_points = prod(grid_base(g) .^ grid_Rs(g))
    if D == 1
        print(io, "DiscretizedGrid{$D}($(grid_base(g)^grid_Rs(g)[1]) points)")
    else
        print(io, "DiscretizedGrid{$D}($(join(grid_base(g) .^ grid_Rs(g), "×")) points)")
    end
end


end
