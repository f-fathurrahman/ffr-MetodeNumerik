# Thomas algorithm for tridigonal system
# a: lower diagonal
# b: main diagonal
# c: upper diagonal
# f: rhs
function linsolve_trid_per!(
    a::Vector{Float64},
    b_::Vector{Float64},
    c::Vector{Float64},
    f::Vector{Float64},
    x::Vector{Float64}
)

    b = copy(b_)
    N = length(b)

    v = zeros(N)
    v[1] = a[1]
    v[N] = c[N]

    # Modified b
    b[1] = b[1] - a[1]
    b[N] = b[N] - c[N]

    Xone = zeros(N)
    Xtwo = zeros(N)

    linsolve_trid!(a, b, c, f, Xone)
    linsolve_trid!(a, b, c, v, Xtwo)

    Xs = ( Xone[1] + Xone[N] )/( 1 + Xtwo[1] + Xtwo[N] )

    for k in 1:N
        x[k] = Xone[k] - Xtwo[k]*Xs
    end

    return
end