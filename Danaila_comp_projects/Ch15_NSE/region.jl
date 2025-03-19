struct Region
    dx::Float64
    dy::Float64
    Lx::Float64
    Ly::Float64
    Nxm::Int64
    Nym::Int64
    idx_ip::Vector{Int64}
    idx_im::Vector{Int64}
    idx_jp::Vector{Int64}
    idx_jm::Vector{Int64}
    idx_ic::Vector{Int64}
    idx_jc::Vector{Int64}
end

function Region(
    Lx::Float64, Ly::Float64,
    Nx::Int64, Ny::Int64
)

    # 2D grid
    Nxm = Nx - 1
    Nym = Ny - 1  # number of cells
    dx = Lx/Nxm
    dy = Ly/Nym

    # Cell indices
    idx_ic = collect(1:Nxm)
    idx_jc = collect(1:Nym)  # indices of cells

    idx_ip = idx_ic .+ 1
    idx_ip[Nxm] = 1

    idx_jp = idx_jc .+ 1
    idx_jp[Nym] = 1  # indices for periodicity

    idx_im = idx_ic .- 1
    idx_im[1] = Nxm

    idx_jm = idx_jc .- 1
    idx_jm[1] = Nym  # cell 1 = cell nxm+1

    return Region(
        dx,
        dy,
        Lx,
        Ly,
        Nxm,
        Nym,
        idx_ip,
        idx_im,
        idx_jp,
        idx_jm,
        idx_ic,
        idx_jc,
    )

end