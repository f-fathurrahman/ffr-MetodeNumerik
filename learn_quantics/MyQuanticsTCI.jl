module MyQuanticsTCI

using MyTensorCrossInterpolation
import MyTensorCrossInterpolation as TCI
import MyQuanticsGrids as QG

import LinearAlgebra: rank
import Base: sum

export quanticscrossinterpolate, evaluate, sum, integral
export cachedata, quanticsfouriermpo

import BitIntegers
import BitIntegers: UInt256, UInt512, UInt1024

BitIntegers.@define_integers 2048 MyInt2048 MyUInt2048

struct QuanticsTensorCI2{ValueType}
    tci::MyTensorCrossInterpolation.TensorCI2{ValueType}
    grid::QG.Grid
    quanticsfunction::TCI.CachedFunction{ValueType}
end

function evaluate(
    qtci::QuanticsTensorCI2{ValueType},
    indices::Union{Array{Int},NTuple{N,Int}}
)::ValueType where {N,ValueType}
    bitlist = QG.grididx_to_quantics(qtci.grid, Tuple(indices))
    return MyTensorCrossInterpolation.evaluate(qtci.tci, bitlist)
end

function evaluate(qtci::QuanticsTensorCI2{V}, indices::Int...)::V where {V}
    return evaluate(qtci, collect(indices)::Vector{Int})
end

function (qtci::QuanticsTensorCI2{V})(indices)::V where {V}
    return evaluate(qtci, indices)
end

function (qtci::QuanticsTensorCI2{V})(indices::Int...)::V where {V}
    return evaluate(qtci, indices...)
end

function sum(qtci::QuanticsTensorCI2{V})::V where {V}
    return sum(qtci.tci)
end

function integral(qtci::QuanticsTensorCI2{V})::V where {V}
    return sum(qtci) * prod(QG.grid_step(qtci.grid))
end


function cachedata(qtci::QuanticsTensorCI2{V}) where {V}
    return Dict(
            QG.quantics_to_origcoord(qtci.grid, k) => v
            for (k, v) in TCI.cachedata(qtci.quanticsfunction)
        )
end

function quanticscrossinterpolate(
    ::Type{ValueType},
    f,
    grid::QG.Grid{n},
    initialpivots::Union{Nothing,AbstractVector{<:AbstractVector}}=nothing;
    nrandominitpivot=5,
    kwargs...
) where {ValueType,n}
    qlocaldimensions = QG.localdimensions(grid)

    qf_ = (n == 1
           ? q -> f(only(QG.quantics_to_origcoord(grid, q)))
           : q -> f(QG.quantics_to_origcoord(grid, q)...))

    maxlinearindex = prod(BigInt.(qlocaldimensions))
    if maxlinearindex < typemax(UInt256)
        keytype = UInt256
    elseif maxlinearindex < typemax(UInt512)
        keytype = UInt512
    elseif maxlinearindex < typemax(UInt1024)
        keytype = UInt1024
    elseif maxlinearindex < typemax(MyUInt2048)
        keytype = MyUInt2048
    else
        keytype = BigInt
    end
    if keytype === BigInt
        @warn "Using BigInt as key type. This will lead to significant memory usage and performance degradation."
    end
    qf = TCI.CachedFunction{ValueType, keytype}(qf_, qlocaldimensions)

    qinitialpivots = (initialpivots === nothing
                      ? [ones(Int, length(qlocaldimensions))]
                      : [QG.grididx_to_quantics(grid, Tuple(p)) for p in initialpivots])

    # For stabity
    kwargs_ = Dict{Symbol,Any}(kwargs)
    if !(:nsearchglobalpivot ∈ keys(kwargs))
        kwargs_[:nsearchglobalpivot] = 5
    end
    if !(:strictlynested ∈ keys(kwargs))
        kwargs_[:strictlynested] = false
    end

    # random initial pivot
    for _ in 1:nrandominitpivot
        pivot = [rand(1:d) for d in qlocaldimensions]
        push!(
            qinitialpivots,
            MyTensorCrossInterpolation.optfirstpivot(qf, qlocaldimensions, pivot)
        )
    end

    qtt, ranks, errors = MyTensorCrossInterpolation.crossinterpolate2(
        ValueType, qf, qlocaldimensions, qinitialpivots; kwargs_...)
    return QuanticsTensorCI2{ValueType}(qtt, grid, qf), ranks, errors
end


function quanticscrossinterpolate(
    ::Type{ValueType},
    f,
    xvals::AbstractVector{<:AbstractVector},
    initialpivots::Union{Nothing,AbstractVector{<:AbstractVector}}=nothing;
    unfoldingscheme::Symbol=:interleaved,
    nrandominitpivot=5,
    kwargs...
) where {ValueType}
    localdimensions = log2.(length.(xvals))
    if !allequal(localdimensions)
        throw(ArgumentError(
            "This method only supports grids with equal number of points in each direction. If you need a different grid, please use index_to_quantics and quantics_to_index and determine the index ordering yourself."))
    elseif !all(isinteger.(localdimensions))
        throw(ArgumentError("This method only supports grid sizes that are powers of 2."))
    end

    n = length(localdimensions)
    R = Int(first(localdimensions))
    grid = QG.DiscretizedGrid{n}(R, Tuple(minimum.(xvals)), Tuple(maximum.(xvals)); unfoldingscheme=unfoldingscheme, includeendpoint=true)

    return quanticscrossinterpolate(ValueType, f, grid, initialpivots; nrandominitpivot=nrandominitpivot, kwargs...)
end

@doc raw"""
    function quanticscrossinterpolate(
        ::Type{ValueType},
        f,
        xvals::AbstractVector,
        initialpivots::AbstractVector=[1];
        kwargs...
    ) where {ValueType}

Interpolate a function ``f(x)`` as a quantics tensor train. This is an overload for 1d functions. For an explanation of arguments and return type, see the documentation of the main overload.
"""
function quanticscrossinterpolate(
    ::Type{ValueType},
    f,
    xvals::AbstractVector,
    initialpivots::AbstractVector=[1];
    nrandominitpivot=5,
    kwargs...
) where {ValueType}
    return quanticscrossinterpolate(
        ValueType,
        f,
        [xvals],
        [initialpivots];
        nrandominitpivot=nrandominitpivot,
        kwargs...)
end

@doc raw"""
    function quanticscrossinterpolate(
        ::Type{ValueType},
        f,
        size::NTuple{d,Int},
        initialpivots::AbstractVector{<:AbstractVector}=[ones(Int, d)];
        unfoldingscheme::Symbol=:interleaved,
        kwargs...
    ) where {ValueType,d}

Interpolate a function ``f(\mathbf{x})`` as a quantics tensor train. This overload automatically constructs a Grid object using the information contained in `size`. Here, the `i`th argument runs from `1` to `size[i]`.
"""
function quanticscrossinterpolate(
    ::Type{ValueType},
    f,
    size::NTuple{d,Int},
    initialpivots::AbstractVector{<:AbstractVector}=[ones(Int, d)];
    unfoldingscheme::Symbol=:interleaved,
    kwargs...
) where {ValueType,d}
    localdimensions = log2.(size)
    if !allequal(localdimensions)
        throw(ArgumentError(
            "This method only supports grids with equal number of points in each direction. If you need a different grid, please use index_to_quantics and quantics_to_index and determine the index ordering yourself."))
    elseif !all(isinteger.(localdimensions))
        throw(ArgumentError("This method only supports grid sizes that are powers of 2."))
    end

    R = Int(first(localdimensions))
    grid = QG.InherentDiscreteGrid{d}(R; unfoldingscheme=unfoldingscheme)
    return quanticscrossinterpolate(ValueType, f, grid, initialpivots; kwargs...)
end

@doc raw"""
    function quanticscrossinterpolate(
        ::Type{ValueType},
        f,
        size::NTuple{d,Int},
        initialpivots::AbstractVector{<:AbstractVector}=[ones(Int, d)];
        unfoldingscheme::Symbol=:interleaved,
        kwargs...
    ) where {ValueType,d}

Interpolate a Tensor ``F`` as a quantics tensor train. For an explanation of arguments, etc., see the documentation of the main overload.
"""
function quanticscrossinterpolate(
    F::Array{ValueType,d},
    initialpivots::AbstractVector{<:AbstractVector}=[ones(Int, d)];
    kwargs...
) where {ValueType,d}
    return quanticscrossinterpolate(
        ValueType,
        (i...) -> F[i...],
        size(F),
        initialpivots;
        kwargs...)
end


module fourierimpl
import MyTensorCrossInterpolation as TCI
struct LagrangePolynomials{T}
    grid::Vector{T}
    baryweights::Vector{T}
end

function (P::LagrangePolynomials{T})(alpha::Int, x::T)::T where {T}
    if abs(x - P.grid[alpha+1]) >= 1e-14
        return prod(x .- P.grid) * P.baryweights[alpha+1] / (x - P.grid[alpha+1])
    else
        return one(T)
    end
end

function getChebyshevGrid(K::Int)::LagrangePolynomials{Float64}
    chebgrid = 0.5 * (1.0 .- cospi.((0:K) / K))
    baryweights = [
        prod(j == m ? 1.0 : 1.0 / (chebgrid[j+1] - chebgrid[m+1]) for m in 0:K)
        for j in 0:K
    ]
    return LagrangePolynomials{Float64}(chebgrid, baryweights)
end

function dftcoretensor(
    P::LagrangePolynomials{Float64},
    alpha::Int, beta::Int, sigma::Int, tau::Int;
    sign::Float64
)::ComplexF64
    x = (sigma + P.grid[beta+1]) / 2
    return P(alpha, x) * cispi(2 * sign * x * tau)
end

end


@doc raw"""
    function quanticsfouriermpo(
        R::Int;
        sign::Float64=-1.0,
        maxbonddim::Int=12,
        tolerance::Float64=1e-14,
        K::Int=25,
        method::Symbol=:SVD,
        normalize::Bool=true
    )::TCI.TensorTrain{ComplexF64}

Generate a quantics Fourier transform operator in tensor train form. When contracted with a quantics tensor train ``F_{\boldsymbol{\sigma}}`` representing a function, the result will be the fourier transform of the function in quantics tensor train form, ``\tilde{F}_{\boldsymbol{\sigma}'} = \sum_{\boldsymbol{\sigma}} F_{\boldsymbol{\sigma}} \exp(-2\pi i (k_{\boldsymbol{\sigma'}}-1) (m_{\boldsymbol{\sigma}} - 1)/M)``, where ``k_{\boldsymbol{\sigma}} = \sum_{\ell=1}^R 2^{R-\ell} \sigma_\ell``, ``m_{\boldsymbol{\sigma}'} = \sum_{\ell=1}^R 2^{R-\ell} \sigma'_\ell``, and ``M=2^R``.

!!! note "Index ordering"
    Before the Fourier transform, the left most index corresponds to ``\sigma_1``, which describes the largest length scale, and the right most index corresponds to ``\sigma_R``, which describes the smallest length scale.
    The indices ``\sigma_1' \ldots \sigma_{R}'`` in the fourier transformed QTT are aligned in the *inverse* order; that is,  the left most index corresponds to ``\sigma'_R``, which describes the smallest length scale.
    This allows construction of an operator with small bond dimension (see reference 1). If necessary, a call to `TCI.reverse(tt)` can restore large-to-small index ordering.

The Fourier transform operator is implemented using a direct analytic construction of the tensor train by Chen and Lindsey (see reference 2). The tensor train thus obtained is then re-compressed to the user-given bond dimension and tolerance.

Arguments:
- `R`: number of bits of the fourier transform.
- `sign`: sign in the exponent ``\exp(2i\pi \times \mathrm{sign} \times (k_{\boldsymbol{\sigma'}}-1) (x_{\boldsymbol{\sigma}}-1)/M)``, usually ``\pm 1``.
- `maxbonddim`: bond dimension to compress the operator to. From observations, `maxbonddim = 12` is generally big enough to reach an accuracy of `1e-12`.
- `tolerance`: tolerance of the TT compression. Note that the error in the fourier transform is generally a bit larger than this error tolerance.
- `K`: bond dimension of the TT before compression, i.e. number of basis functions to approximate the Fourier transform with (see reference 2). The TT will become inaccurate for `K < 22`; higher values may be necessary for very high precision.
- `method`: method with which to compress the TT. Choose between `:SVD` and `:CI`.
- `normalize`: whether or not to normalize the operator as an isometry.

!!! details "References"
    1. [J. Chen, E. M. Stoudenmire, and S. R. White, Quantum Fourier Transform Has Small Entanglement, PRX Quantum 4, 040318 (2023).](https://link.aps.org/doi/10.1103/PRXQuantum.4.040318)
    2. [J. Chen and M. Lindsey, Direct Interpolative Construction of the Discrete Fourier Transform as a Matrix Product Operator, arXiv:2404.03182.](http://arxiv.org/abs/2404.03182)
"""
function quanticsfouriermpo(
    R::Int;
    sign::Float64=-1.0,
    tolerance::Float64=1e-14,
    maxbonddim::Int=12,
    K::Int=25, method::Symbol=:SVD,
    normalize::Bool=true
)::TCI.TensorTrain{ComplexF64}
    P = fourierimpl.getChebyshevGrid(K)
    A = [
        fourierimpl.dftcoretensor(P, alpha, beta, sigma, tau; sign)
        for alpha in 0:K, tau in [0, 1], sigma in [0, 1], beta in 0:K
    ]
    Afirst = reshape(sum(A, dims=1), (1, 2, 2, K + 1))
    Alast = reshape(A[:, :, :, 1], (K + 1, 2, 2, 1))
    tt = TCI.TensorTrain{ComplexF64,4}([Afirst, fill(A, R - 2)..., Alast])
    TCI.compress!(tt, method; tolerance, maxbonddim)
    if normalize
        for t in tt.sitetensors
            t ./= sqrt(2.0)
        end
    end
    return tt
end

function swaphalves!(tt::TCI.TensorTrain{V,3}) where {V}
    tt.sitetensors[1][:, :, :] = tt.sitetensors[1][:, [2, 1], :]
    nothing
end

function swaphalves(tt::TCI.AbstractTensorTrain{V}) where {V}
    ttcopy = TCI.tensortrain(deepcopy(TCI.sitetensors(tt)))
    swaphalves!(ttcopy)
    return ttcopy
end


end
