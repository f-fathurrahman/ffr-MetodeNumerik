# Thomas algorithm for tridigonal system
# a: lower diagonal
# b: main diagonal
# c: upper diagonal
# f: rhs
function linsolve_trid!(
    a::Vector{Float64},
    b::Vector{Float64},
    c::Vector{Float64},
    f::Vector{Float64},
    x::Vector{Float64}
)

    N = length(b)

    β = zeros(N)
    β[1] = b[1]
    for k in 2:N
        β[k] = b[k] - c[k-1]/β[k-1]*a[k]
    end

    γ = zeros(N)    
    γ[1] = f[1]/β[1]
    for k in 2:N
        γ[k] = ( f[k] - a[k]*γ[k-1] )/β[k]
    end
    println("γ = ", γ)

    x[N] = γ[N]
    for k in (N-1):-1:1
        x[k] = γ[k] - c[k]/β[k]*x[k+1]
        println("k = ", k)
        println("x[k] = ", x[k])
    end

    return
end